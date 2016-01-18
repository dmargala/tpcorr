#!/usr/bin/env python

import argparse
import numpy as np
import scipy.interpolate
import h5py

import astropy.time
import astropy.coordinates
import astropy.units as u
import astropy.units.imperial
import astropy.units.cds

import bossdata
import specsim

from pointing import *
from guider import *
from acceptance_model import *

# helper function for shifting angles to desired range
def normalize_angle(angle):
    while angle <= -180:
        angle += 360
    while angle > 180:
        angle -= 360
    return angle

# calculate throughputs
def calculate_target_offsets(plate, mjd, guide_wlen=5400*u.Angstrom, std_wlen=5400.*u.Angstrom,
                             offset_wlen=4000*u.Angstrom, steps_per_exposure=5, wlen_grid_steps=15, save_guide_plots=False):

    print 'Calculating corrections for {} observed on MJD {}'.format(plate, mjd)
    
    finder = bossdata.path.Finder()
    mirror = bossdata.remote.Manager()

    # Get the list of exposures used in this observation's coadd from a spec lite file.
    spec_name = finder.get_spec_path(plate, mjd, fiber=1, lite=True)
    spec_file = bossdata.spec.SpecFile(mirror.get(spec_name))

    # Read the first b1 raw science exposure to find this plate's plug map.
    raw = spec_file.get_raw_image(0, 'blue', finder=finder, mirror=mirror)
    plug_map = raw.read_plug_map()

    # Look up the plate design pointing from the raw header and convert to
    # an index A,B,C,... -> 0,1,2,...
    pointing_label = raw.header['POINTING'].strip()
    pointing_index = ord(pointing_label) - ord('A')
    
    # Initialize a pointing object for this plate's sky location.
    ra_center = float(plug_map['raCen']) * u.deg
    dec_center = float(plug_map['decCen']) * u.deg
    print 'Plate center is RA={:.3f}, DEC={:.3f} for {}-{}'.format(ra_center, dec_center, plate, pointing_label)
    pointing = Pointing(ra_center, dec_center)
    
    # Find the nominal observing temperature and time that this plate's holes are drilled for.
    # apo = specsim.transform.observatories['APO']
    apo = astropy.coordinates.EarthLocation.of_site('apo')

    design_temp = float(plug_map['temp'])*u.deg_C
    design_pressure = None # Calculate based on elevation and temperature
    design_ha = float(plug_map['ha'].split()[pointing_index]) * u.deg
    midnight = astropy.time.Time(mjd, format='mjd', scale='tai', location=apo)
    design_time = specsim.transform.adjust_time_to_hour_angle(midnight, ra_center, design_ha)
    design_tai = design_time.mjd * 86400.
    print 'Holes drilled for T={:.1f} and HA={:.1f} (TAI={:.1f})'.format(design_temp, design_ha, design_tai)
    
    # Find this plate's guide stars.
    plugging = plug_map['PLUGMAPOBJ']
    guide_fibers = plugging['holeType'] == 'GUIDE'

    guide_ra, guide_dec = plugging['ra'][guide_fibers], plugging['dec'][guide_fibers]
    guide_targets = astropy.coordinates.ICRS(guide_ra * u.deg, guide_dec * u.deg)
    
    # Calculate the nominal guide fiber positions.
    guide_x0, guide_y0, _, _ = pointing.transform(guide_targets, design_tai, guide_wlen,
                                                  design_temp, design_pressure)
    
    # Find this plate's offset fibers. We have to use spAll for this since the plug map does
    # not record the design wavelengths.
    # spAll = bossdata.meta.Database(lite=False)
    # offset_fibers = spAll.select_all(
    #     where='PLATE={} and MJD={} and LAMBDA_EFF={}'
    #     .format(plate, mjd, offset_wlen.to(u.Angstrom).value),
    #     what='FIBER,OBJTYPE,PLUG_RA,PLUG_DEC,XFOCAL,YFOCAL')

    # this is a little hacky for now, spAll isn't available at this point in the pipeline so
    # this is a first pass at removing the dependency by using data in spFrame files
    b1_frame_name = finder.get_plate_path(plate, spec_file.get_exposure_name(0, 'blue', 'spFrame'))
    b1_frame = bossdata.plate.FrameFile(mirror.get(b1_frame_name))
    # load frame file for second spectograph
    spec2_name = finder.get_spec_path(plate, mjd, fiber=501, lite=True)
    spec2_file = bossdata.spec.SpecFile(mirror.get(spec2_name))
    b2_frame_name = finder.get_plate_path(plate, spec2_file.get_exposure_name(0, 'blue', 'spFrame'))
    b2_frame = bossdata.plate.FrameFile(mirror.get(b2_frame_name))
    # combine plugmap info
    frame_plugmap = astropy.table.vstack([b1_frame.plug_map, b2_frame.plug_map])
    # Find offset fibers
    offset_fibers_mask = frame_plugmap['LAMBDA_EFF'] == offset_wlen.to(u.Angstrom).value
    offset_fibers = frame_plugmap[offset_fibers_mask]

    offset_xfocal = offset_fibers['XFOCAL'] * u.mm
    offset_yfocal = offset_fibers['YFOCAL'] * u.mm
    offset_targets = astropy.coordinates.ICRS(
        ra=offset_fibers['RA'] * u.deg,
        dec=offset_fibers['DEC'] * u.deg)
    print 'Plate has {:d} guide fibers and {:d} offset targets.'.format(
        len(guide_targets), np.count_nonzero(offset_targets))

    # Calculate the nominal science fiber positions. These will not match XFOCAL, YFOCAL
    # exactly since we do not exactly replicate the IDL transforms, but they should be
    # close (within ~0.2 arcsec) and we only use offsets calculated consistently with
    # transform() in the following.
    offset_x0, offset_y0, offset_alt, offset_az = pointing.transform(
        offset_targets, design_tai, offset_wlen, design_temp, design_pressure)
    
    # Calculate where the offset target fibers would have been positioned if they were
    # designed for the same wavelength as the standard stars.
    offset_x0_std, offset_y0_std, _, _ = pointing.transform(
        offset_targets, design_tai, std_wlen, design_temp, design_pressure)
    
    # Initialize the wavelength grid to use for calculating corrections.
    wlen_grid = np.linspace(3500., 10500., wlen_grid_steps)[:, np.newaxis] * u.Angstrom
    
    # Initialize guided target centroid list
    guided_centroids = []

    # Initialize seeing array
    obs_seeing = []
    obs_ha = []
    obs_pressure = []
    obs_temperature = []

    # Precompute the conversion from inches of Hg to kPa.
    pconv = (1 * u.cds.mmHg * u.imperial.inch / u.mm).to(u.kPa).value
    
    # Loop over exposures
    for exp_index in range(spec_file.num_exposures):

        # Open the b1 frame for this exposure, to access its metadata.
        b1_frame_name = finder.get_plate_path(
            plate, spec_file.get_exposure_name(exp_index, 'blue', 'spFrame'))
        b1_frame = bossdata.plate.FrameFile(mirror.get(b1_frame_name))
        exp_id = b1_frame.exposure_id

        # Lookup this exposure's observing time, seeing, and temperature.
        tai_beg, tai_end = b1_frame.header['TAI-BEG'], b1_frame.header['TAI-END']
        tai_mid = 0.5 * (tai_beg + tai_end)
        
        ha = normalize_angle(pointing.hour_angle(tai_mid).to(u.deg).value)*u.deg
        
        seeing = b1_frame.header['SEEING50'] * u.arcsec
        try:
            temperature = b1_frame.header['AIRTEMP'] * u.deg_C
        except ValueError,e:
            temperature = design_temp
        try:
            pressure = b1_frame.header['PRESSURE'] * pconv * u.kPa
        except ValueError,e:
            pressure = design_pressure

        print 'Exp[{:02d}] #{:08d} seeing {:.3f}, T={:+5.1f}, P={:.1f}, TAI {:.1f} ({:+7.3f} days, HA {:+.1f})'.format(
            exp_index, exp_id, seeing, temperature, pressure, tai_mid, (tai_mid - design_tai)/86400., ha.to(u.deg))
        
        obs_seeing.append(seeing)
        obs_ha.append(ha)
        obs_pressure.append(pressure)
        obs_temperature.append(temperature)
        
        # Create time steps covering this exposure.
        tai_steps = np.linspace(tai_beg, tai_end, steps_per_exposure)

        # Calculate the actual guide target positions on the focal plane without any guiding.
        guide_x, guide_y, _, _ = pointing.transform(
            guide_targets[:, np.newaxis], tai_steps, guide_wlen, temperature, pressure)
        
        # Solve for the optimal guider corrections.
        guider = Guider(guide_x0, guide_y0, guide_x, guide_y)
        if save_guide_plots:
            guider.plot(tai_steps, field_radius=340 * u.mm, zoom=5000.,
                        fiber_radius=0.1 * u.arcsec * pointing.platescale,
                        save='guide-{}-{}-{}.png'.format(plate, mjd, exp_index))

        # Calculate the offset target paths on the focal plane without any guiding, for the
        # actual observing conditions.
        offset_x, offset_y, _, _ = pointing.transform(
            offset_targets[:, np.newaxis, np.newaxis], tai_steps, wlen_grid, temperature, pressure)
        
        # Apply guiding corrections to estimate the actual offset target paths during the exposure.
        guided_centroids.append(guider.correct(offset_x, offset_y))

    results = {
        'fibers':offset_fibers['FIBERID'], 'num_exposures':spec_file.num_exposures,
        'obs_seeing':obs_seeing, 'obs_ha':obs_ha, 'obs_temperature':obs_temperature, 'obs_pressure':obs_pressure, 
        'wlen_grid':wlen_grid, 'steps_per_exposure':steps_per_exposure,
        'guided_centroids':guided_centroids,
        'offset_x0':offset_x0, 'offset_y0':offset_y0,
        'offset_x0_std':offset_x0_std, 'offset_y0_std':offset_y0_std
        }

    return results

def calculate_corrections(offsets, seeing_wlen=5400.*u.Angstrom, platescale=217.7358*u.mm/u.deg):

    # Initialize telescope model
    sdss_25m = Telescope(diameter=2.5*u.m, obscuration_area_fraction=0.27, plate_scale=platescale)

    # Extract precomputed centroid positions
    offset_x0 = offsets['offset_x0']
    offset_y0 = offsets['offset_y0']
    offset_x0_std = offsets['offset_x0_std']
    offset_y0_std = offsets['offset_y0_std']
    guided_centroids = offsets['guided_centroids']

    corrections = np.empty(
        (offsets['num_exposures'], len(offsets['fibers']), len(offsets['wlen_grid']), offsets['steps_per_exposure']),
        dtype=float)

    obs_seeing = offsets['obs_seeing']

    for exp_index in range(offsets['num_exposures']):

        seeing = obs_seeing[exp_index]

        psf = sdss_25m.get_atmospheric_psf(seeing_wlen, seeing, gauss=False)
        fiber_diameter = 2*u.arcsec
        acceptance_model = sdss_25m.calculate_fiber_acceptance(fiber_diameter, psf)

        guided_x, guided_y = guided_centroids[exp_index]

        # Calculate centroid offsets for each offset target, relative to its nominal fiber center.
        offset = np.sqrt(
            (guided_x - offset_x0[:, np.newaxis, np.newaxis])**2 +
            (guided_y - offset_y0[:, np.newaxis, np.newaxis])**2)
        
        # Calculate centroid offsets for each offset target, relative to where its fiber center would
        # be if it were designed for the same wavelength as the standard stars.
        offset_std = np.sqrt(
            (guided_x - offset_x0_std[:, np.newaxis, np.newaxis])**2 +
            (guided_y - offset_y0_std[:, np.newaxis, np.newaxis])**2)

        # Calculate the acceptance fractions for both sets of centroid offsets.
        acceptance = acceptance_model(0.5*(offset / platescale).to(u.arcsec).value)
        acceptance_std = acceptance_model(0.5*(offset_std / platescale).to(u.arcsec).value)
        
        # Calculate the acceptance fraction ratios, tabulated for each offset target, wavelength and time.
        # The ratio calculated this way gives the correction of eqn (13).
        corrections[exp_index] = acceptance_std / acceptance

    # Average the correction over each exposure time slice.
    avg_corrections = np.mean(np.mean(corrections, axis=-1), axis=0)

    return corrections, avg_corrections

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--plate', type=int, default=6641,
        help='plate id')
    parser.add_argument('-m', '--mjd', type=str, default=None,
        help='observation mjd')
    parser.add_argument('--output', type=str, default=None,
        help='output filename')
    parser.add_argument('--wlen-grid-steps', type=int, default=15,
        help='Number of wlen grid steps to use for correction calculation (between 3500-10500 incl.)')
    parser.add_argument('--save-guide-plots', action='store_true',
        help='Save per exposure guide plots.')
    parser.add_argument('--save-debug-data', action='store_true',
        help='Save per exposure offset data.')
    args = parser.parse_args()

    finder = bossdata.path.Finder()
    mirror = bossdata.remote.Manager()

    # If no mjd is provided, try to guess. For plates with multiple pluggings, print list of mjds and exit.
    if args.mjd is None:
        mjd_list = bossdata.meta.get_plate_mjd_list(
            args.plate, finder=finder, mirror=mirror)
        if len(mjd_list) == 0:
            print('Plate {0:d} has not been observed with good data quality.'.format(
                args.plate))
            return -1
        elif len(mjd_list) > 1:
            print('Plate {0:d} has been observed on MJDs {1:s}.'.format(
                args.plate, ','.join(map(str, mjd_list))))
            print('Select one of these using the --mjd command-line argument.')
            return -1
        else:
            args.mjd = mjd_list[0]

    # open output file
    filename = 'corrections-%s-%s.hdf5' % (str(args.plate), str(args.mjd)) if args.output is None else args.output
    outfile = h5py.File(filename, 'w')

    # Calculate offsets for individual exposures
    offsets = calculate_target_offsets(args.plate, args.mjd, wlen_grid_steps=args.wlen_grid_steps,
        save_guide_plots=args.save_guide_plots)

    # Calculate average correction from individual exposure offsets
    corrections, avg_corrections = calculate_corrections(offsets)

    # Save corrections to output file
    outfile.create_dataset('wave', data=offsets['wlen_grid'])
    outfile.create_group(str(args.plate))
    outfile[str(args.plate)].create_group(str(args.mjd))

    for fiber, avg_correction in zip(offsets['fibers'], avg_corrections):
        dset = outfile.create_dataset('/'.join(map(str,[args.plate,args.mjd,fiber])), data=avg_correction, dtype='f4')
    outfile.close()

    if args.save_debug_data:
        filename = 'debug_' + filename
        outfile = h5py.File(filename, 'w')
        outfile.create_dataset('wave', data=offsets['wlen_grid'])
        plate_grp = outfile.create_group(str(args.plate))
        platemjd_grp = plate_grp.create_group(str(args.mjd))

        exp_corrections = np.mean(corrections, axis=-1)

        platemjd_grp.create_dataset('offset_x0', data=offsets['offset_x0'].to(u.mm).value)
        platemjd_grp.create_dataset('offset_y0', data=offsets['offset_y0'].to(u.mm).value)
        platemjd_grp.create_dataset('offset_x0_std', data=offsets['offset_x0_std'].to(u.mm).value)
        platemjd_grp.create_dataset('offset_y0_std', data=offsets['offset_y0_std'].to(u.mm).value)
        platemjd_grp.create_dataset('fibers', data=offsets['fibers'])
        platemjd_grp.create_dataset('obs_seeing', data=[seeing.to(u.arcsec).value for seeing in offsets['obs_seeing']])
        platemjd_grp.create_dataset('obs_ha', data=[ha.to(u.degree).value for ha in offsets['obs_ha']])
        platemjd_grp.create_dataset('obs_pressure', data=[pressure.to(u.kPa).value for pressure in offsets['obs_pressure']])
        platemjd_grp.create_dataset('obs_temperature', data=[temperature.to(u.deg_C).value for temperature in offsets['obs_temperature']])

        for exp_index in range(offsets['num_exposures']):
            exp_grp = platemjd_grp.create_group(str(exp_index))
            guided_x, guided_y = offsets['guided_centroids'][exp_index]
            exp_guided_x = np.mean(guided_x, axis=-1)
            exp_guided_y = np.mean(guided_y, axis=-1)

            for i,fiber in enumerate(offsets['fibers']):
                fiber_grp = exp_grp.create_group(str(fiber))

                fiber_grp.attrs['x0'] = offsets['offset_x0'][i].to(u.mm).value
                fiber_grp.attrs['y0'] = offsets['offset_y0'][i].to(u.mm).value
                fiber_grp.attrs['x0_std'] = offsets['offset_x0_std'][i].to(u.mm).value
                fiber_grp.attrs['y0_std'] = offsets['offset_y0_std'][i].to(u.mm).value
                fiber_grp.create_dataset('correction', data=exp_corrections[exp_index,i])
                fiber_grp.create_dataset('guided_x', data=exp_guided_x[i].to(u.mm).value)
                fiber_grp.create_dataset('guided_y', data=exp_guided_y[i].to(u.mm).value)


if __name__ == "__main__":
    main()

