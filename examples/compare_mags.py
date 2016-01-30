#!/usr/bin/env python
"""
"""
import argparse
import os

import h5py
import numpy as np

import astropy.table
from astropy.io import fits

import scipy.interpolate

import bossdata

import specsim

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tpcorr', type=str, default=None,
        help='throughput correction filename, required')
    parser.add_argument('--output', type=str, default=None,
        help='output filename')
    parser.add_argument('--blue-path', type=str, default='/sas/dr12/boss',
        help='path to blue reduction')
    parser.add_argument('--blue-version', type=str, default='test',
        help='blue reduction version')
    args = parser.parse_args()

    # Open connection spAll db
    meta_db = bossdata.meta.Database(lite=False, verbose=True)

    validiation_plates = (6130,6131,6135,6136,6147,6155,6157,6290,6293,6296,6297,6298,6307,6506,6509,6590,6681,6734,6816,6986)

    # Summarize various target sample selections used in paper
    sample_names = ['Failed quasars', 'Spec. standards', 'Offset standards', 
        'Failed quasars 2', 'Offset standards 2']
    sample_selections = [
        'LAMBDA_EFF=4000 and OBJTYPE="QSO" and CLASS="STAR"',
        'LAMBDA_EFF=5400 and OBJTYPE="SPECTROPHOTO_STD" and CLASS="STAR"',
        'LAMBDA_EFF=4000 and CLASS="STAR" and ANCILLARY_TARGET2=(1<<20)',
        'LAMBDA_EFF=4000 and OBJTYPE="QSO" and CLASS="STAR"',
        'LAMBDA_EFF=4000 and OBJTYPE="SPECTROPHOTO_STD" and CLASS="STAR"',
    ]

    # Require targets were on "validation" plates
    validiation_plates_str = ','.join(['{}'.format(plate) for plate in validiation_plates])
    valid_selection = 'ZWARNING=0 and PLATE in ({})'.format(validiation_plates_str)

    bad_chunks = ('boss35','boss36','boss37','boss38')
    bad_chunks_str = ','.join(['"{}"'.format(chunk) for chunk in bad_chunks])
    valid2_selection = 'ZWARNING=0 and CHUNK in ({})'.format(bad_chunks_str)
    
    # Loop over target samples
    what = 'PLATE,MJD,FIBER,PSFMAG_1,PSFMAG_2,PSFMAG_3'
    sql_prefix = 'SELECT {} FROM meta'.format(what)
    target_lists = []
    for sample_name, sample_selection in zip(sample_names, sample_selections):
        sample_nums = {}
        # Count the number of targets in this sample+category
        if '2' in sample_name:
            sql = sql_prefix + ' WHERE {} and {}'.format(sample_selection, valid2_selection)
        else:
            sql = sql_prefix + ' WHERE {} and {}'.format(sample_selection, valid_selection)
        meta_db.cursor.execute(sql)
        rows = meta_db.cursor.fetchall()
        print len(rows)
        target_lists.append(rows)

    # Open the throughput correction file
    tpcorr = h5py.File(args.tpcorr, 'r')
    tpcorr_wave = tpcorr['wave'].value

    # Initialize finders
    dr12_finder = bossdata.path.Finder(sas_path='/sas/dr12/boss', redux_version='v5_7_0')
    blue_finder = bossdata.path.Finder(sas_path=args.blue_path, redux_version=args.blue_version)
    mirror = bossdata.remote.Manager()

    bands = ['g','r','i']

    # Need to correct sdss imaging magnitudes 
    # More more details see: http://www.sdss3.org/dr8/algorithms/fluxcal.php#SDSStoAB
    ab_minus_sdss = {'g':0.036, 'r':0.015, 'i':0.013}
    # ab_minus_sdss = {'g':0.012, 'r':0.010, 'i':0.028}

    final_sample_names = [sample_names[1], sample_names[2], 'Corr. offset standards', 
        'Spec. offset standards', sample_names[0], 'Corr. failed quasars', sample_names[3], sample_names[4]]
    final_target_list_indices = [1,2,2,2,0,0,3,4]

    for name, target_list_index in zip(final_sample_names, final_target_list_indices):

        # Allocate arrays for storing magnitude/color (g-r, r-i, g, r, i) calculations
        imaging_mags = np.ma.empty((len(target_lists[target_list_index]), 5))
        syn_mags = np.ma.empty_like(imaging_mags)
        imaging_mags[:] = np.ma.masked
        syn_mags[:] = np.ma.masked

        # Loop over targets
        for i,target in enumerate(target_lists[target_list_index]):
            # Parse target row tuple
            plate,mjd,fiber,psfmag_1,psfmag_2,psfmag_3 = target
            # Assign imaging magnitudes to dictionary
            sdss_mags = {'g':psfmag_1, 'r':psfmag_2, 'i':psfmag_3}
            # Correct sdss imaging magnitudes
            corrected_sdss_mags = {band: ab_minus_sdss[band] + sdss_mags[band] for band in bands}

            # Get this target's spectrum
            # For the "spectroscopic" offset standard sample, use the special "blue" reduction
            if name == 'Spec. offset standards':
                blue_plate_filename = blue_finder.get_plate_spec_path(plate, mjd)
                blue_platefile = bossdata.plate.PlateFile(mirror.get(blue_plate_filename))
                fibers = np.array([fiber], dtype=int)
                data = blue_platefile.get_valid_data(fibers, fiducial_grid=True, use_ivar=True)[0]
            # Otherwise use the dr12 lite spec file
            else:
                dr12_spec_name = dr12_finder.get_spec_path(plate, mjd, fiber=fiber, lite=True)
                dr12_spec = bossdata.spec.SpecFile(mirror.get(dr12_spec_name))
                data = dr12_spec.get_valid_data(fiducial_grid=True, use_ivar=True)
            wlen,flux,ivar = data['wavelength'][:],data['flux'][:],data['ivar'][:]

            # For the "corrected" samples, look up and play the target's throughut correction
            if name[:4] == 'Corr':
                # Load this target's correction
                correction = tpcorr['{}/{}/{}'.format(plate,mjd,fiber)].value
                # Create an interpolated correction function
                correction_interp = scipy.interpolate.interp1d(tpcorr_wave, correction, kind='linear', 
                    bounds_error=False, fill_value=np.ma.masked)
                # Sample the interpolated correction using the observation's wavelength grid
                resampled_correction = correction_interp(wlen)
                # Apply the correction to the observed flux and ivar
                flux *= resampled_correction
                ivar /= (resampled_correction**2)

            # Calculate synthetic ab magnitudes for this spectrum
            spectrum = specsim.spectrum.SpectralFluxDensity(wlen, flux)
            syn_mag = spectrum.getABMagnitudes()

            # If there's a problem with any of the magnitudes, throw the target away
            missing_syn_mag = False
            for band in bands:
                if syn_mag[band] == None:
                    missing_syn_mag = True
            if missing_syn_mag:
                continue

            # Save imaging and synthetic magnitudes
            imaging_mags[i,2] = corrected_sdss_mags['g']
            imaging_mags[i,3] = corrected_sdss_mags['r']
            imaging_mags[i,4] = corrected_sdss_mags['i']

            syn_mags[i,2] = syn_mag['g']
            syn_mags[i,3] = syn_mag['r']
            syn_mags[i,4] = syn_mag['i']

        # Calculate colors
        imaging_mags[:,0] = imaging_mags[:,2] - imaging_mags[:,3]
        imaging_mags[:,1] = imaging_mags[:,3] - imaging_mags[:,4]
        syn_mags[:,0] = syn_mags[:,2] - syn_mags[:,3]
        syn_mags[:,1] = syn_mags[:,3] - syn_mags[:,4]

        # Compare imaging and synthetic magnitudes
        delta = syn_mags - imaging_mags
        delta_mean = np.ma.mean(delta, axis=0)

        # Normalized Median Absolute Deviation: a "Robust" version of standard deviation
        # https://en.wikipedia.org/wiki/Median_absolute_deviation 
        delta_median = np.ma.median(delta, axis=0)
        delta_mad = np.ma.median(np.ma.abs(delta - delta_median), axis=0)
        delta_nmad = 1.4826*delta_mad
        
        # Print summary table
        summary = '   '.join(['{: .3f} {: .3f}'.format(mean, disp) for mean, disp in zip(delta_mean, delta_nmad)])
        print '{:<25}: {}'.format(name, summary)

if __name__ == '__main__':
    main()
