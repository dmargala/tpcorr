#!/usr/bin/env python
"""
"""
import argparse
import os

import h5py
import numpy as np

from astropy.io import fits

import scipy.interpolate
import scipy.stats.mstats as mstats
import scipy.signal

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.size': 14})
mpl.rcParams.update({'savefig.dpi': 100})
mpl.rcParams.update({'savefig.bbox': 'tight'})
import matplotlib.pyplot as plt

import bossdata

def plot_median_flux_ratio(wlen, red_ratio, blue_ratio, red_label=None, blue_label=None, show_quantiles=False):

    # Calculate median and 1-sigma quantiles
    quantiles = [scipy.stats.norm.cdf(sigma) for sigma in [-1, 0, 1]]
    red_quantiles = [mstats.mquantiles(red_ratio, q, axis=0)[0] for q in quantiles]
    blue_quantiles = [mstats.mquantiles(blue_ratio, q, axis=0)[0] for q in quantiles]

    # Plot median
    plt.plot(wlen, red_quantiles[1], color='red', label=red_label)
    plt.plot(wlen, blue_quantiles[1], color='blue', label=blue_label)

    if show_quantiles:
        plt.fill_between(wlen, red_quantiles[0], red_quantiles[-1], color='red', alpha=0.3, lw=0)
        plt.fill_between(wlen, blue_quantiles[0], blue_quantiles[-1], color='blue', alpha=0.3, lw=0)
        
    plt.xlim(wlen[0], wlen[-1])
    plt.ylabel('Median Flux Ratio')
    plt.xlabel(r'Observed Wavelength $\lambda$ $(\mathrm{\AA})$')
    plt.grid()

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tpcorr', type=str, default=None,
        help='throughput correction filename, required')
    parser.add_argument('--output', type=str, default=None,
        help='output filename')
    parser.add_argument('--validation', action='store_true',
        help='validation reduction')
    args = parser.parse_args()

    # Open the throughput correction file
    tpcorr = h5py.File(args.tpcorr, 'r')
    tpcorr_wave = np.squeeze(tpcorr['wave'].value)

    sdss_finder = bossdata.path.Finder(sas_path='/sas/dr12/sdss', redux_version='26')
    sdss_finder_103 = bossdata.path.Finder(sas_path='/sas/dr12/sdss', redux_version='103')
    sdss_finder_104 = bossdata.path.Finder(sas_path='/sas/dr12/sdss', redux_version='104')
    dr12_finder = bossdata.path.Finder(sas_path='/sas/dr12/boss', redux_version='v5_7_0')
    blue_finder = bossdata.path.Finder(sas_path='/sas/dr12/boss', redux_version='test')

    mirror = bossdata.remote.Manager()

    def get_sdss_spec(plate, mjd, fiber):
        try:
            sdss_spec_name = sdss_finder.get_spec_path(plate, mjd, fiber=fiber, lite=True)
            return bossdata.spec.SpecFile(mirror.get(sdss_spec_name))
        except RuntimeError:
            try:
                sdss_spec_name = sdss_finder_103.get_spec_path(plate, mjd, fiber=fiber, lite=True)
                return bossdata.spec.SpecFile(mirror.get(sdss_spec_name))
            except RuntimeError:
                try:
                    sdss_spec_name = sdss_finder_104.get_spec_path(plate, mjd, fiber=fiber, lite=True)
                    return bossdata.spec.SpecFile(mirror.get(sdss_spec_name))
                except RuntimeError:
                    print 'Could not find sdss spectrum ({}-{}-{}) in redux 26, 203, or 104.'.format(plate,mjd,fiber)
        return None

    if args.validation:
        # Open connection spAll db
        meta_db = bossdata.meta.Database(lite=False, verbose=True)
        # Failed quasars
        sample_selection = 'LAMBDA_EFF=4000 and OBJTYPE="QSO" and CLASS="STAR"'
        # Require targets were on "validation" plates
        validiation_plates = (6130,6131,6135,6136,6147,6155,6157,6290,6293,6296,6297,6298,6307,6506,6509,6590,6681,6734,6816,6986)
        validiation_plates_str = ','.join(['{}'.format(plate) for plate in validiation_plates])
        valid_selection = 'ZWARNING=0 and PLATE in ({})'.format(validiation_plates_str)
        what = 'PLATE,MJD,FIBER'
        sort = 'PLATE,MJD,FIBER'
        where = '{} and {}'.format(sample_selection, valid_selection)
        targets = meta_db.select_all(what, where, sort, max_rows=0)
    else:
        # Find repeat quasar observations using DR12Q (55752 cutoff is approx DR9 cutoff to compare with Paris et al composite)
        dr12q_db = bossdata.meta.Database(quasar_catalog=True, quasar_catalog_name='DR12Q')
        what = 'PLATE,MJD,FIBER,Z_VI,PLATE_DR7,MJD_DR7,FIBERID_DR7'
        where = 'SDSS_DR7!=0 and BAL_FLAG_VI=0 and ZWARNING&(1<<7)=0 and MJD <= 55752'
        sort = 'PLATE,MJD,FIBER'
        targets = dr12q_db.select_all(what, where, sort, max_rows=0)

    num_targets = len(targets)

    print num_targets

    # Use the "fiducial" wavelength grid
    fiducial_loglam = bossdata.spec.fiducial_loglam
    fiducial_wavelength = 10**fiducial_loglam
    num_pixels = len(fiducial_loglam)

    # Allocate masked arrays to store fluxes and ivars
    ref_composite_flux = np.ma.empty((num_targets, num_pixels))
    ref_composite_ivar = np.ma.empty_like(ref_composite_flux)
    ref_composite_flux[:] = np.ma.masked
    ref_composite_ivar[:] = np.ma.masked

    dr12_composite_flux = np.ma.empty_like(ref_composite_flux)
    dr12_composite_ivar = np.ma.empty_like(ref_composite_flux)
    dr12_composite_flux[:] = np.ma.masked
    dr12_composite_ivar[:] = np.ma.masked

    tpcorr_composite_flux = np.ma.empty_like(ref_composite_flux)
    tpcorr_composite_ivar = np.ma.empty_like(ref_composite_flux)
    tpcorr_composite_flux[:] = np.ma.masked
    tpcorr_composite_ivar[:] = np.ma.masked

    num_good_spectra = 0

    for i,target in enumerate(targets):
        plate, mjd, fiber = target['PLATE'],target['MJD'],target['FIBER']
        if i and i % 100 == 0:
            print i, plate, mjd, fiber, num_good_spectra

        # Get the DR12 spectrum
        dr12_spec_name = dr12_finder.get_spec_path(plate, mjd, fiber=fiber, lite=True)
        dr12_spec = bossdata.spec.SpecFile(mirror.get(dr12_spec_name, timeout=None, progress_min_size=0))
        
        # Make sure this is an offset target
        lambda_eff = dr12_spec.hdulist[2].read_column('LAMBDA_EFF')[0]
        if lambda_eff != 4000:
            # print '{}-{}-{}: LAMBDA_EFF!=4000'.format(plate, mjd, fiber)
            continue

        # Check that correction exists for this target
        tpcorr_key = '/'.join(map(str,[plate,mjd,fiber]))
        if tpcorr_key not in tpcorr:
            print '{}-{}-{}: No correction in tpcorr file.'.format(plate, mjd, fiber)
            continue

        # Get reference fiducial spectrum
        if args.validation:
            spplate_filename = blue_finder.get_plate_spec_path(plate, mjd)
            platefile = bossdata.plate.PlateFile(mirror.get(spplate_filename))
            fibers = np.array([fiber], dtype=int)
            data = platefile.get_valid_data(fibers, fiducial_grid=True, use_ivar=True)[0]
        else:
            sdss_plate, sdss_mjd, sdss_fiberid = target['PLATE_DR7'],target['MJD_DR7'],target['FIBERID_DR7']
            sdss_spec = get_sdss_spec(sdss_plate, sdss_mjd, sdss_fiberid)

            if sdss_spec is None:
                continue
            data = sdss_spec.get_valid_data(fiducial_grid=True, use_ivar=True)

        # All checks passed
        num_good_spectra += 1

        wlen,flux,ivar = data['wavelength'][:],data['flux'][:],data['ivar'][:]
        ref_composite_flux[i] = flux
        ref_composite_ivar[i] = ivar

        # Save flux and ivar
        data = dr12_spec.get_valid_data(fiducial_grid=True, use_ivar=True)
        wlen,flux,ivar = data['wavelength'][:],data['flux'][:],data['ivar'][:]
        dr12_composite_flux[i] = flux
        dr12_composite_ivar[i] = ivar

        # Read the target's throughput correction vector
        correction = tpcorr[tpcorr_key].value
        # Create an interpolated correction function
        correction_interp = scipy.interpolate.interp1d(tpcorr_wave, correction, kind='linear', 
            bounds_error=False, fill_value=np.ma.masked)
        # Sample the interpolated correction using the observation's wavelength grid
        resampled_correction = correction_interp(wlen)
        # Apply the correction to the observed flux and ivar
        tpcorr_composite_flux[i] = flux * resampled_correction
        tpcorr_composite_ivar[i] = ivar / (resampled_correction**2)

    print num_targets, num_good_spectra

    fig = plt.figure(figsize=(16,4))
    dr12_ratio = dr12_composite_flux/ref_composite_flux
    tpcorr_ratio = tpcorr_composite_flux/ref_composite_flux
    plot_median_flux_ratio(fiducial_wavelength, dr12_ratio, tpcorr_ratio, show_quantiles=True)

    if args.validation:
        plt.ylim(0.7, 1.3)
        plt.savefig(args.output)
    else:
        plt.ylim(0, 2)
        plt.savefig(args.output)

    # if args.output:
    #     # save target list with sn column
    #     outfile = h5py.File(args.output, 'w')

    #     outfile.create_dataset('sdss_flux', data=ref_composite_flux)
    #     outfile.create_dataset('sdss_ivar', data=ref_composite_ivar)
    #     outfile.create_dataset('sdss_mask', data=ref_composite_flux.mask)
    #     outfile.create_dataset('dr12_flux', data=dr12_composite_flux)
    #     outfile.create_dataset('dr12_ivar', data=dr12_composite_ivar)
    #     outfile.create_dataset('dr12_mask', data=dr12_composite_flux.mask)
    #     outfile.create_dataset('tpcorr_flux', data=tpcorr_composite_flux)
    #     outfile.create_dataset('tpcorr_ivar', data=tpcorr_composite_ivar)
    #     outfile.create_dataset('tpcorr_mask', data=tpcorr_composite_flux.mask)
    #     outfile.create_dataset('wavelength', data=fiducial_wavelength)
    #     outfile.attrs['ntargets'] = num_targets

    #     outfile.close()


if __name__ == '__main__':
    main()
