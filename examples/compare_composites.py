#!/usr/bin/env python
"""
"""
import argparse
import os

import h5py
import numpy as np

from astropy.io import fits
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt

import bossdata


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tpcorr', type=str, default=None,
        help='throughput correction filename, required')
    parser.add_argument('--target-list', type=str, default=None,
        help='target id string (plate-mjd-fiberid), must specify --spec-dir')
    parser.add_argument('--output', type=str, default=None,
        help='output filename')
    parser.add_argument('--sdss-lookup', type=str, default=None,
        help='sdss lookup table')
    parser.add_argument('--validation', action='store_true',
        help='validation reduction')
    parser.add_argument('--use-bossdata', action='store_true',
        help='use bossdata to load specfiles')
    args = parser.parse_args()

    # Open the throughput correction file
    tpcorr = h5py.File(args.tpcorr, 'r')
    tpcorr_wave = tpcorr['wave'].value

    target_ids = []

    sdss_finder = bossdata.path.Finder(sas_path='/sas/dr12/sdss', redux_version='26')
    dr12_finder = bossdata.path.Finder(sas_path='/sas/dr12/boss', redux_version='v5_7_0')
    mirror = bossdata.remote.Manager()

    sdss_lookup = {}
    if not args.validation:
        with open(args.sdss_lookup) as infile:
            for line in infile:
                fields = line.split()
                sdss_lookup[fields[0]] = '%s-%s-%s' % (fields[1], fields[2], fields[3])

    with open(args.target_list) as infile:
        for line in infile:
            boss_target_id = line.split()[0]

            if not args.validation:
                try:
                    sdss_target_id = sdss_lookup[boss_target_id]
                except KeyError:
                    continue

            try:
                tpcorr_key = '/'.join(boss_target_id.split('-'))
                correction = tpcorr[tpcorr_key].value
            except KeyError:
                pass

            target_ids.append(boss_target_id)

    print len(target_ids)

    min_fid_index = 0
    max_fid_index = 4800
    npixels = max_fid_index - min_fid_index

    fiducial_loglam = bossdata.spec.fiducial_loglam
    num_pixels = len(fiducial_loglam)

    # norm_min_index = np.round(get_fiducial_pixel_index_offset(np.log10(args.norm_min))).astype(int)-min_fid_index
    # norm_max_index = np.round(get_fiducial_pixel_index_offset(np.log10(args.norm_max))).astype(int)-min_fid_index

    sdss_composite_flux = np.ma.empty((len(target_ids), num_pixels))
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

    dr12_norms = np.zeros(len(target_ids))
    tpcorr_norms = np.zeros(len(target_ids))
    sdss_norms = np.zeros(len(target_ids))

    plate_ids = np.empty(len(target_ids))
    fiber_ids = np.empty(len(target_ids))

    nboss = 0
    ntpcorr = 0
    nsdss = 0
    counter = 0

    sdss_finder = bossdata.path.Finder(sas_path='/sas/dr12/sdss', redux_version='26')
    dr12_finder = bossdata.path.Finder(sas_path='/sas/dr12/boss', redux_version='v5_7_0')
    blue_finder = bossdata.path.Finder(sas_path='/sas/dr12/boss', redux_version='test')

    for i, target_id in enumerate(target_ids):

        if i and i % 50 == 0:
            print i, target_id

        plate, mjd, fiberid = [int(field) for field in target_id.split('-')]

        plate_ids[i] = plate
        fiber_ids[i] = mjd

        if args.use_bossdata:
            dr12_spec_name = dr12_finder.get_spec_path(plate, mjd, fiber=fiberid, lite=True)
            dr12_spec = bossdata.spec.SpecFile(mirror.get(dr12_spec_name, timeout=None))
            data = dr12_spec.get_valid_data(fiducial_grid=True, use_ivar=True)
            wlen,flux,ivar = data['wavelength'][:],data['flux'][:],data['ivar'][:]

            dr12_composite_flux[i] = flux
            dr12_composite_ivar[i] = ivar
            dr12_norms[i] = 1

            # Read the target's throughput correction vector
            tpcorr_key = '%s/%s/%s' % (plate, mjd, fiberid)
            correction = tpcorr[tpcorr_key].value

            # Create an interpolated correction function
            correction_interp = interp1d(tpcorr_wave, correction, kind='linear', 
                bounds_error=False, fill_value=np.ma.masked)

            # Sample the interpolated correction using the observation's wavelength grid
            resampled_correction = correction_interp(wlen)

            # Apply the correction to the observed flux and ivar
            tpcorr_composite_flux[i] = flux*resampled_correction
            tpcorr_composite_ivar[i] = ivar / (resampled_correction**2)
            tpcorr_norms[i] = 1

            # Repeat for sdss observation
            sdss_target_id = sdss_lookup[target_id]
            sdss_plate, sdss_mjd, sdss_fiberid = [int(field) for field in sdss_target_id.split('-')]
            sdss_spec_name = sdss_finder.get_spec_path(sdss_plate, sdss_mjd, fiber=sdss_fiberid, lite=True)
            sdss_spec = bossdata.spec.SpecFile(mirror.get(sdss_spec_name, timeout=None))
            data = sdss_spec.get_valid_data(fiducial_grid=True, use_ivar=True)
            wlen,flux,ivar = data['wavelength'][:],data['flux'][:],data['ivar'][:]

            sdss_composite_flux[i] = flux
            sdss_composite_ivar[i] = ivar
            sdss_norms[i] = 1

        else:
            dr12_spec_name = dr12_finder.get_spec_path(plate, mjd, fiber=fiberid, lite=True)
            # Load the target's combined spectrum
            try:
                spec = fits.open(mirror.get(dr12_spec_name, timeout=None))
            except IOError:
                raise IOError('Error opening spec file: %s' % dr12_spec_name)

            flux = spec[1].data.field('flux')
            ivar = spec[1].data.field('ivar')
            and_mask = spec[1].data.field('and_mask')
            ivar[and_mask > 0] = 0

            # z = spec[2].data.field('z')[0]
            loglam = spec[1].data.field('loglam')
            wavelength = np.power(10, loglam)
            spec.close()

            pixel_offsets = bossdata.spec.get_fiducial_pixel_index(wavelength)
            pixel_offset_indices = np.round(pixel_offsets).astype(int)
            composite_indices = pixel_offset_indices - min_fid_index

            valid_pixels = (composite_indices < npixels) & (composite_indices >= 0)

            pixel_offset = bossdata.spec.get_fiducial_pixel_index(wavelength[0])

            dr12_composite_flux[i, composite_indices[valid_pixels]] = flux[valid_pixels]
            dr12_composite_ivar[i, composite_indices[valid_pixels]] = ivar[valid_pixels]
            dr12_norms[i] = 1

            # Read the target's throughput correction vector
            tpcorr_key = '%s/%s/%s' % (plate, mjd, fiberid)
            correction = tpcorr[tpcorr_key].value

            # Create an interpolated correction function
            correction_interp = interp1d(np.squeeze(tpcorr_wave), correction, kind='linear', 
                bounds_error=False, fill_value=np.ma.masked)

            # Sample the interpolated correction using the observation's wavelength grid
            resampled_correction = correction_interp(wavelength)

            # Apply the correction to the observed flux and ivar
            tpcorr_composite_flux[i, composite_indices[valid_pixels]] = (flux*resampled_correction)[valid_pixels]
            tpcorr_composite_ivar[i, composite_indices[valid_pixels]] = (ivar / (resampled_correction**2))[valid_pixels]
            tpcorr_norms[i] = 1

            if args.validation:

                spplate_filename = blue_finder.get_plate_spec_path(plate, mjd)
                # spplate_filename = os.path.join(args.boss_dir, version, '%04d' % plate, 'spPlate-%04d-%5d.fits' % (plate, mjd))
                # Load the target's combined spectrum
                try:
                    spplate = fits.open(mirror.get(spplate_filename, timeout=None))
                except IOError:
                    raise IOError('Error opening spPlate file: %s' % spplate_filename)

                flux = spplate[0].data[fiberid-1]
                ivar = spplate[1].data[fiberid-1]
                and_mask = spplate[2].data[fiberid-1]

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
                sdss_norms[i] = 1

            else:

                # Repeat for sdss observation
                sdss_target_id = sdss_lookup[target_id]
                sdss_plate, sdss_mjd, sdss_fiberid = [int(field) for field in sdss_target_id.split('-')]
                sdss_spec_name = sdss_finder.get_spec_path(sdss_plate, sdss_mjd, fiber=sdss_fiberid, lite=True)

                # Load the target's combined spectrum
                try:
                    spec = fits.open(mirror.get(sdss_spec_name, timeout=None))
                except IOError:
                    raise IOError('Error opening spec file: %s' % sdss_spec_name)

                sdss_spec = bossdata.spec.SpecFile(mirror.get(sdss_spec_name, timeout=None))
                
                flux = spec[1].data.field('flux')
                ivar = spec[1].data.field('ivar')
                and_mask = spec[1].data.field('and_mask')
                ivar[and_mask > 0] = 0

                # z = spec[2].data.field('z')[0]
                loglam = spec[1].data.field('loglam')
                wavelength = np.power(10, loglam)
                spec.close()

                pixel_offsets = bossdata.spec.get_fiducial_pixel_index(wavelength)
                pixel_offset_indices = np.round(pixel_offsets).astype(int)
                composite_indices = pixel_offset_indices - min_fid_index

                valid_pixels = (composite_indices < npixels) & (composite_indices >= 0)

                pixel_offset = bossdata.spec.get_fiducial_pixel_index(wavelength[0])

                sdss_composite_flux[i, composite_indices[valid_pixels]] = flux[valid_pixels]
                sdss_composite_ivar[i, composite_indices[valid_pixels]] = ivar[valid_pixels]
                sdss_norms[i] = 1

    print nboss, ntpcorr, nsdss

    print sdss_composite_flux.shape

    if args.output:
        # save target list with sn column
        outfile = h5py.File(args.output, 'w')

        outfile.create_dataset('sdss_flux', data=sdss_composite_flux)
        outfile.create_dataset('sdss_ivar', data=sdss_composite_ivar)
        outfile.create_dataset('sdss_mask', data=sdss_composite_flux.mask)
        outfile.create_dataset('sdss_norm', data=sdss_norms)
        outfile.create_dataset('dr12_flux', data=dr12_composite_flux)
        outfile.create_dataset('dr12_ivar', data=dr12_composite_ivar)
        outfile.create_dataset('dr12_mask', data=dr12_composite_flux.mask)
        outfile.create_dataset('dr12_norm', data=dr12_norms)
        outfile.create_dataset('tpcorr_flux', data=tpcorr_composite_flux)
        outfile.create_dataset('tpcorr_ivar', data=tpcorr_composite_ivar)
        outfile.create_dataset('tpcorr_mask', data=tpcorr_composite_flux.mask)
        outfile.create_dataset('tpcorr_norm', data=tpcorr_norms)
        outfile.create_dataset('plate_ids', data=plate_ids)
        outfile.create_dataset('fiber_ids', data=fiber_ids)
        outfile.create_dataset('wavelength', data=10**fiducial_loglam)
        outfile.attrs['ntargets'] = len(target_ids)

        outfile.close()


if __name__ == '__main__':
    main()
