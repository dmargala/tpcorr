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

import matplotlib.pyplot as plt

import bossdata

def ratios(repeats, ymin=0.5, ymax=1.5, xlabel=r'Observed Wavelength $\lambda$ $(\AA)$', nfilter=None, show_quantiles=False):

    wave_min = repeats['wavelength'].value[0]
    wave_max = repeats['wavelength'].value[-1]    
    wave = repeats['wavelength'].value
    
    boss_flux = np.ma.masked_array(repeats['dr12_flux'].value, mask=repeats['dr12_mask'].value)
    sdss_flux = np.ma.masked_array(repeats['sdss_flux'].value, mask=repeats['sdss_mask'].value)
    tpcorr_flux = np.ma.masked_array(repeats['tpcorr_flux'].value, mask=repeats['tpcorr_mask'].value)
            
    boss_ratio = (boss_flux)/sdss_flux
    tpcorr_ratio = (tpcorr_flux)/sdss_flux
    
    if nfilter:
        boss_ratio = scipy.signal.medfilt(boss_ratio, [1, nfilter])
        tpcorr_ratio = scipy.signal.medfilt(tpcorr_ratio, [1, nfilter])
    
    print wave.shape, boss_flux.shape

    quantiles = [scipy.stats.norm.cdf(sigma) for sigma in [-1, 0, 1]]
    boss_quantiles = [mstats.mquantiles(boss_ratio, q, axis=0)[0] for q in quantiles]
    tpcorr_quantiles = [mstats.mquantiles(tpcorr_ratio, q, axis=0)[0] for q in quantiles]
    print boss_quantiles[0].shape, boss_quantiles[-1].shape

    plt.plot(wave, boss_quantiles[1], color='red', label='BOSS')
    plt.plot(wave, tpcorr_quantiles[1], color='blue', label='Corrected BOSS')

    if show_quantiles:
        plt.fill_between(wave, boss_quantiles[0], boss_quantiles[-1], color='red', alpha=.3, lw=0)
        plt.fill_between(wave, tpcorr_quantiles[0], tpcorr_quantiles[-1], color='blue', alpha=.3, lw=0)
        
        
    plt.ylim(ymin,ymax)
    plt.xlim(wave_min, wave_max)
    plt.ylabel('Median Flux Ratio')
    plt.xlabel(xlabel)
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
    tpcorr_wave = tpcorr['wave'].value

    sdss_finder = bossdata.path.Finder(sas_path='/sas/dr12/sdss', redux_version='26')
    sdss_finder_103 = bossdata.path.Finder(sas_path='/sas/dr12/sdss', redux_version='103')
    sdss_finder_104 = bossdata.path.Finder(sas_path='/sas/dr12/sdss', redux_version='104')
    dr12_finder = bossdata.path.Finder(sas_path='/sas/dr12/boss', redux_version='v5_7_0')
    blue_finder = bossdata.path.Finder(sas_path='/sas/dr12/boss', redux_version='test')

    mirror = bossdata.remote.Manager()

    # Find repeat quasar observations using DR12Q
    dr12q_db = bossdata.meta.Database(quasar_catalog=True, quasar_catalog_name='DR12Q')
    what = 'PLATE,MJD,FIBER,Z_VI,PLATE_DR7,MJD_DR7,FIBERID_DR7'
    where = 'SDSS_DR7!=0 and BAL_FLAG_VI=0 and MJD <= 55752'
    sort = 'PLATE,MJD,FIBER'
    repeat_targets = dr12q_db.select_all(what, where, sort, max_rows=0)

    num_repeat_targets = len(repeat_targets)

    repeat4000 = []
    repeat5400 = []

    # bossquery --where 'LAMBDA_EFF=4000 and OBJTYPE="SPECTROPHOTO_STD"' --full --print --verbose --max-rows 0 --save spec_std_4000.txt
    # cat spec_std_4000.txt | tail -n +2 | cut -d' ' -f 1 | uniq
    blue_std_plates = (7033, 7034, 7035, 7036, 7037, 7038, 7039, 7040, 7041, 7042, 
        7043, 7044, 7045, 7046, 7047, 7048, 7049, 7050, 7051, 7052, 7053, 7054, 7055,
        7056, 7057, 7058, 7059, 7146, 7147, 7148, 7149, 7150, 7151, 7152, 7155, 7159,
        7160, 7161, 7162, 7163, 7164, 7165, 7166, 7167, 7168, 7169, 7183, 7184)


    fiducial_loglam = bossdata.spec.fiducial_loglam
    num_pixels = len(fiducial_loglam)

    sdss_composite_flux = np.ma.empty((num_repeat_targets, num_pixels))
    sdss_composite_ivar = np.ma.empty_like(sdss_composite_flux)
    sdss_composite_flux[:] = np.ma.masked
    sdss_composite_ivar[:] = np.ma.masked

    dr12_composite_flux = np.ma.empty_like(sdss_composite_flux)
    dr12_composite_ivar = np.ma.empty_like(sdss_composite_flux)
    dr12_composite_flux[:] = np.ma.masked
    dr12_composite_ivar[:] = np.ma.masked

    tpcorr_composite_flux = np.ma.empty_like(sdss_composite_flux)
    tpcorr_composite_ivar = np.ma.empty_like(sdss_composite_flux)
    tpcorr_composite_flux[:] = np.ma.masked
    tpcorr_composite_ivar[:] = np.ma.masked

    for i,target in enumerate(repeat_targets):
        plate, mjd, fiber = target['PLATE'],target['MJD'],target['FIBER']
        if i and i % 100 == 0:
            print i, plate, mjd, fiber

        if plate in blue_std_plates:
            continue

        tpcorr_key = '/'.join(map(str,[plate,mjd,fiber]))
        try:
            correction = tpcorr[tpcorr_key].value
            lambda_eff = 4000
        except KeyError:
            lambda_eff = 5400

        dr12_spec_name = dr12_finder.get_spec_path(plate, mjd, fiber=fiber, lite=True)
        dr12_spec = bossdata.spec.SpecFile(mirror.get(dr12_spec_name, timeout=None, progress_min_size=0))

        if dr12_spec.hdulist[2].read_column('LAMBDA_EFF')[0] != lambda_eff:
            print plate, mjd, fiber, lambda_eff, dr12_spec.hdulist[2].read_column('LAMBDA_EFF')[0]

        data = dr12_spec.get_valid_data(fiducial_grid=True, use_ivar=True)
        wlen,flux,ivar = data['wavelength'][:],data['flux'][:],data['ivar'][:]

        dr12_composite_flux[i] = flux
        dr12_composite_ivar[i] = ivar

        # Read the target's throughput correction vector
        if lambda_eff == 4000:
            correction = tpcorr[tpcorr_key].value

            # Create an interpolated correction function
            correction_interp = scipy.interpolate.interp1d(tpcorr_wave, correction, kind='linear', 
                bounds_error=False, fill_value=np.ma.masked)

            # Sample the interpolated correction using the observation's wavelength grid
            resampled_correction = correction_interp(wlen)

            # Apply the correction to the observed flux and ivar
            tpcorr_composite_flux[i] = flux*resampled_correction
            tpcorr_composite_ivar[i] = ivar / (resampled_correction**2)

        if not args.validation:
            # Repeat for sdss observation
            sdss_plate, sdss_mjd, sdss_fiberid = target['PLATE_DR7'],target['MJD_DR7'],target['FIBERID_DR7']
            try:
                sdss_spec_name = sdss_finder.get_spec_path(sdss_plate, sdss_mjd, fiber=sdss_fiberid, lite=True)
                sdss_spec = bossdata.spec.SpecFile(mirror.get(sdss_spec_name, timeout=None, progress_min_size=0))
            except RuntimeError:
                try:
                    sdss_spec_name = sdss_finder_103.get_spec_path(sdss_plate, sdss_mjd, fiber=sdss_fiberid, lite=True)
                    sdss_spec = bossdata.spec.SpecFile(mirror.get(sdss_spec_name, timeout=None, progress_min_size=0))
                except RuntimeError:
                    try:
                        sdss_spec_name = sdss_finder_104.get_spec_path(sdss_plate, sdss_mjd, fiber=sdss_fiberid, lite=True)
                        sdss_spec = bossdata.spec.SpecFile(mirror.get(sdss_spec_name, timeout=None, progress_min_size=0))
                    except RuntimeError:
                        continue

            data = sdss_spec.get_valid_data(fiducial_grid=True, use_ivar=True)
            wlen,flux,ivar = data['wavelength'][:],data['flux'][:],data['ivar'][:]

            sdss_composite_flux[i] = flux
            sdss_composite_ivar[i] = ivar

        else:
            spplate_filename = blue_finder.get_plate_spec_path(plate, mjd)
            # spplate_filename = os.path.join(args.boss_dir, version, '%04d' % plate, 'spPlate-%04d-%5d.fits' % (plate, mjd))
            # Load the target's combined spectrum
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

            spplate.close()

            pixel_offsets = bossdata.spec.get_fiducial_pixel_index(10**loglam)
            pixel_offset_indices = np.round(pixel_offsets).astype(int)
            composite_indices = pixel_offset_indices - min_fid_index

            valid_pixels = (composite_indices < npixels) & (composite_indices >= 0)

            pixel_offset = bossdata.spec.get_fiducial_pixel_index(wavelength[0])

            sdss_composite_flux[i, composite_indices[valid_pixels]] = flux[valid_pixels]
            sdss_composite_ivar[i, composite_indices[valid_pixels]] = ivar[valid_pixels]

    print sdss_composite_flux.shape

    if args.output:
        # save target list with sn column
        outfile = h5py.File(args.output, 'w')

        outfile.create_dataset('sdss_flux', data=sdss_composite_flux)
        outfile.create_dataset('sdss_ivar', data=sdss_composite_ivar)
        outfile.create_dataset('sdss_mask', data=sdss_composite_flux.mask)
        outfile.create_dataset('dr12_flux', data=dr12_composite_flux)
        outfile.create_dataset('dr12_ivar', data=dr12_composite_ivar)
        outfile.create_dataset('dr12_mask', data=dr12_composite_flux.mask)
        outfile.create_dataset('tpcorr_flux', data=tpcorr_composite_flux)
        outfile.create_dataset('tpcorr_ivar', data=tpcorr_composite_ivar)
        outfile.create_dataset('tpcorr_mask', data=tpcorr_composite_flux.mask)
        outfile.create_dataset('wavelength', data=10**fiducial_loglam)
        outfile.attrs['ntargets'] = num_repeat_targets

        outfile.close()


if __name__ == '__main__':
    main()
