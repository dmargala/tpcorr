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

def mad(arr):
    """
    Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    # arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.ma.median(arr, axis=0)
    return np.ma.median(np.ma.abs(arr - med), axis=0)

def nmad(arr):
    return 1.4826*mad(arr)

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tpcorr', type=str, default=None,
        help='throughput correction filename, required')
    parser.add_argument('--output', type=str, default=None,
        help='output filename')
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

    dr12_finder = bossdata.path.Finder(sas_path='/sas/dr12/boss', redux_version='v5_7_0')
    blue_finder = bossdata.path.Finder(sas_path='/sas/dr12/boss', redux_version='test')
    mirror = bossdata.remote.Manager()


    bands = ['g','r','i']
    ab_minus_sdss = {'g':0.036, 'r':0.015, 'i':0.013}
    # ab_minus_sdss = {'g':0.012, 'r':0.010, 'i':0.028}

    final_sample_names = [sample_names[1], sample_names[2], 'Corr. offset standards', 
        'Spec. offset standards', sample_names[0], 'Corr. failed quasars', sample_names[3], sample_names[4]]
    final_target_list_indices = [1,2,2,2,0,0,3,4]

    for name, target_list_index in zip(final_sample_names, final_target_list_indices):

        imaging_grid = np.empty((len(target_lists[target_list_index]),5))
        syn_grid = np.empty_like(imaging_grid)

        for i,target in enumerate(target_lists[target_list_index]):

            plate,mjd,fiber,psf1,psf2,psf3 = target

            sdss_mags = {'g':psf1, 'r':psf2, 'i':psf3}
            corrected_sdss_mags = {band: ab_minus_sdss[band] + sdss_mags[band] for band in bands}

            if name == 'Spec. offset standards':
                spplate_filename = blue_finder.get_plate_spec_path(plate, mjd)
                try:
                    spplate = fits.open(mirror.get(spplate_filename, timeout=None))
                except IOError:
                    raise IOError('Error opening spPlate file: %s' % spplate_filename)

                flux = spplate[0].data[fiber-1]
                ivar = spplate[1].data[fiber-1]
                and_mask = spplate[2].data[fiber-1]

                coeff0 = spplate[0].header['COEFF0']
                coeff1 = spplate[0].header['COEFF1']
                loglam = coeff0 + coeff1*np.arange(0, len(flux))
                wlen = 10**loglam

            else:
                dr12_spec_name = dr12_finder.get_spec_path(plate, mjd, fiber=fiber, lite=True)
                dr12_spec = bossdata.spec.SpecFile(mirror.get(dr12_spec_name, progress_min_size=0))
                data = dr12_spec.get_valid_data(fiducial_grid=True, use_ivar=True)
                wlen,flux,ivar = data['wavelength'][:],data['flux'][:],data['ivar'][:]

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

            spectrum = specsim.spectrum.SpectralFluxDensity(wlen, flux)
            syn_mag = spectrum.getABMagnitudes()

            imaging_grid[i,0] = (corrected_sdss_mags['g']-corrected_sdss_mags['r'])
            imaging_grid[i,1] = (corrected_sdss_mags['r']-corrected_sdss_mags['i'])
            imaging_grid[i,2] = corrected_sdss_mags['g']
            imaging_grid[i,3] = corrected_sdss_mags['r']
            imaging_grid[i,4] = corrected_sdss_mags['i']

            syn_grid[i,0] = (syn_mag['g']-syn_mag['r'])
            syn_grid[i,1] = (syn_mag['r']-syn_mag['i'])
            syn_grid[i,2] = syn_mag['g']
            syn_grid[i,3] = syn_mag['r']
            syn_grid[i,4] = syn_mag['i']

        delta_mean = np.mean(syn_grid - imaging_grid, axis=0)
        delta_disp = nmad(syn_grid - imaging_grid)
        
        summary = '   '.join(['{: .3f} {: .3f}'.format(mean, disp) for mean, disp in zip(delta_mean, delta_disp)])

        print '{:<25}: {}'.format(name, summary)

if __name__ == '__main__':
    main()
