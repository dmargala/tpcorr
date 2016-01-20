## Observation correction model

import numpy as np

import astropy.time
import astropy.coordinates
import astropy.units as u
import astropy.units.imperial
import astropy.units.cds

import bossdata
import specsim

import tpcorr.pointing
import tpcorr.guider
import tpcorr.acceptance_model

class Observation(object):
    def __init__(self, plate, mjd, guide_wlen=5400 * u.Angstrom, offset_wlen=4000 * u.Angstrom,
        std_wlen=5400 * u.Angstrom, wlen_grid_steps=15, steps_per_exposure=5):
        print 'Calculating corrections for {} observed on MJD {}'.format(plate, mjd)

        self.plate = plate
        self.mjd = mjd
        self.guide_wlen = guide_wlen
        self.std_wlen = std_wlen
        self.offset_wlen = offset_wlen
        self.steps_per_exposure = steps_per_exposure
        self.wlen_grid_steps = wlen_grid_steps
        
        self.finder = bossdata.path.Finder()
        self.mirror = bossdata.remote.Manager()

        # Get the list of exposures used in this observation's coadd from a spec lite file.
        spec_name = self.finder.get_spec_path(plate, mjd, fiber=1, lite=True)
        self.spec_file = bossdata.spec.SpecFile(self.mirror.get(spec_name))

        # Read the first b1 raw science exposure to find this plate's plug map.
        raw = self.spec_file.get_raw_image(0, 'blue', finder=self.finder, mirror=self.mirror)
        plug_map = raw.read_plug_map()

        # Look up the plate design pointing from the raw header and convert to
        # an index A,B,C,... -> 0,1,2,...
        pointing_label = raw.header['POINTING'].strip()
        pointing_index = ord(pointing_label) - ord('A')
        
        # Initialize a pointing object for this plate's sky location.
        ra_center = float(plug_map['raCen']) * u.deg
        dec_center = float(plug_map['decCen']) * u.deg
        print 'Plate center is RA={:.3f}, DEC={:.3f} for {}-{}'.format(ra_center, dec_center, plate, pointing_label)
        self.pointing = tpcorr.pointing.Pointing(ra_center, dec_center)
        
        # Find the nominal observing temperature and time that this plate's holes are drilled for.
        design_temp = float(plug_map['temp'])*u.deg_C
        design_pressure = None # Calculate based on elevation and temperature
        design_ha = float(plug_map['ha'].split()[pointing_index]) * u.deg
        midnight = astropy.time.Time(mjd, format='mjd', scale='tai', location=self.pointing.where)
        design_time = specsim.transform.adjust_time_to_hour_angle(midnight, ra_center, design_ha)
        self.design_tai = design_time.mjd * 86400.
        print 'Holes drilled for T={:.1f} and HA={:.1f} (TAI={:.1f})'.format(design_temp, design_ha, self.design_tai)
        
        # Find this plate's guide stars.
        plugging = plug_map['PLUGMAPOBJ']
        guide_fibers = plugging['holeType'] == 'GUIDE'

        guide_ra, guide_dec = plugging['ra'][guide_fibers], plugging['dec'][guide_fibers]
        self.guide_targets = astropy.coordinates.ICRS(guide_ra * u.deg, guide_dec * u.deg)
        
        # Calculate the nominal guide fiber positions.
        self.guide_x0, self.guide_y0, _, _ = self.pointing.transform(
            self.guide_targets, self.design_tai, guide_wlen, design_temp, design_pressure)
        
        # Find this plate's offset fibers. We have to use spAll for this since the plug map does
        # not record the design wavelengths.
        self.plugmap = self.get_plugmap_from_spframes()

        offset_fibers_mask = self.plugmap['LAMBDA_EFF'] == offset_wlen.to(u.Angstrom).value
        offset_fibers = self.plugmap[offset_fibers_mask]

        offset_xfocal = offset_fibers['XFOCAL'] * u.mm
        offset_yfocal = offset_fibers['YFOCAL'] * u.mm
        self.fiber_ids = offset_fibers['FIBERID']
        self.offset_targets = astropy.coordinates.ICRS(ra=offset_fibers['RA'] * u.deg, dec=offset_fibers['DEC'] * u.deg)
        self.num_offset_targets = np.count_nonzero(self.offset_targets)
        print 'Plate has {:d} guide fibers and {:d} offset targets.'.format(len(self.guide_targets), self.num_offset_targets)

        # Calculate the nominal science fiber positions. These will not match XFOCAL, YFOCAL
        # exactly since we do not exactly replicate the IDL transforms, but they should be
        # close (within ~0.2 arcsec) and we only use offsets calculated consistently with
        # transform() in the following.
        self.offset_x0, self.offset_y0, offset_alt, offset_az = self.pointing.transform(
            self.offset_targets, self.design_tai, offset_wlen, design_temp, design_pressure)
        
        # Calculate where the offset target fibers would have been positioned if they were
        # designed for the same wavelength as the standard stars.
        self.offset_x0_std, self.offset_y0_std, _, _ = self.pointing.transform(
            self.offset_targets, self.design_tai, std_wlen, design_temp, design_pressure)
        
        # Initialize the wavelength grid to use for calculating corrections.
        self.wlen_grid = np.linspace(3500., 10500., wlen_grid_steps)[:, np.newaxis] * u.Angstrom
        
        # Initialize guided target centroid list
        self.guided_centroids = []

        # Initialize exposure meta data lists
        self.seeing = np.empty((self.spec_file.num_exposures)) * u.arcsec
        self.ha = np.empty((self.spec_file.num_exposures)) * u.degree
        self.pressure = np.empty((self.spec_file.num_exposures)) * u.kPa
        self.temperature = np.empty((self.spec_file.num_exposures)) * u.deg_C
        self.tai_beg = np.empty((self.spec_file.num_exposures)) # seconds
        self.tai_end = np.empty((self.spec_file.num_exposures)) # seconds

        self.init_exposure_meta(temperature0=design_temp)


    def init_exposure_meta(self, seeing0=1.49 * u.arcsec, pressure0=79.3 * u.kPa, temperature0=5 * u.deg_C):
        # Precompute the conversion from inches of Hg to kPa.
        pconv = (1 * u.cds.mmHg * u.imperial.inch / u.mm).to(u.kPa).value

        # Loop over exposures
        for exp_index in range(self.spec_file.num_exposures):

            # Open the b1 frame for this exposure, to access its metadata.
            b1_frame_name = self.finder.get_plate_path(
                self.plate, self.spec_file.get_exposure_name(exp_index, 'blue', 'spFrame'))
            b1_frame = bossdata.plate.FrameFile(self.mirror.get(b1_frame_name))
            exp_id = b1_frame.exposure_id

            # Lookup this exposure's observing time, seeing, and temperature.
            self.tai_beg[exp_index] = b1_frame.header['TAI-BEG']
            self.tai_end[exp_index] = b1_frame.header['TAI-END']
            tai_mid = 0.5 * (self.tai_beg[exp_index] + self.tai_end[exp_index])
            
            # Convert tai to hour angle
            self.ha[exp_index] = tpcorr.pointing.normalize_angle(
                self.pointing.hour_angle(tai_mid).to(u.deg).value)*u.deg
            
            if b1_frame.header['SEEING50'] == 0:
                print 'Warning: SEEING50=0, using nominal value.'
                self.seeing[exp_index] = seeing0
            else:
                self.seeing[exp_index] = b1_frame.header['SEEING50'] * u.arcsec

            try:
                self.temperature[exp_index] = b1_frame.header['AIRTEMP'] * u.deg_C
            except ValueError, e:
                print 'Warning: AIRTEMP not found, using nominal value:', e
                self.temperature[exp_index] = temperature0
            try:
                self.pressure[exp_index] = b1_frame.header['PRESSURE'] * pconv * u.kPa
            except ValueError, e:
                print 'Warning: PRESSURE not found, using nominal value.', e
                self.pressure[exp_index] = pressure0

            print 'Exp[{:02d}] #{:08d} seeing {:.3f}, T={:+5.1f}, P={:.1f}, TAI {:.1f} ({:+7.3f} days, HA {:+.1f})'.format(
                exp_index, exp_id, self.seeing[exp_index], self.temperature[exp_index], self.pressure[exp_index], 
                tai_mid, (tai_mid - self.design_tai)/86400., self.ha[exp_index])

    def get_plugmap_from_spframes(self):
        # Read frame files for both spectrographs
        frames = []
        for fiber in (1,501):
            spec_name = self.finder.get_spec_path(self.plate, self.mjd, fiber=fiber, lite=True)
            spec_file = bossdata.spec.SpecFile(self.mirror.get(spec_name))
            frame_name = self.finder.get_plate_path(self.plate, spec_file.get_exposure_name(0, 'blue', 'spFrame'))
            frames.append(bossdata.plate.FrameFile(self.mirror.get(frame_name)))
        # Stack frame plugmaps
        return astropy.table.vstack([frame.plug_map for frame in frames])

    def get_exp_centroids(self, exp_index, guide_plot_name=None):
        # Create time steps covering this exposure.
        tai_steps = np.linspace(self.tai_beg[exp_index], self.tai_end[exp_index], self.steps_per_exposure)

        # Calculate the actual guide target positions on the focal plane without any guiding.
        guide_x, guide_y, _, _ = self.pointing.transform(
            self.guide_targets[:, np.newaxis], tai_steps, self.guide_wlen, self.temperature[exp_index], self.pressure[exp_index])
        
        # Solve for the optimal guider corrections.
        guider = tpcorr.guider.Guider(self.guide_x0, self.guide_y0, guide_x, guide_y)
        if guide_plot_name:
            guider.plot(tai_steps, field_radius=340 * u.mm, zoom=5000., 
                fiber_radius=0.1 * u.arcsec * self.pointing.platescale, save=guide_plot_name)

        # Calculate the offset target paths on the focal plane without any guiding, for the actual observing conditions.
        offset_x, offset_y, _, _ = self.pointing.transform(
            self.offset_targets[:, np.newaxis, np.newaxis], tai_steps, self.wlen_grid, self.temperature[exp_index], self.pressure[exp_index])
        
        # Apply guiding corrections to estimate the actual offset target paths during the exposure.
        return guider.correct(offset_x, offset_y)


    def get_corrections(self, seeing_wlen=5400.*u.Angstrom):

        # Precompute wlen ratio for wavelength dependent seeing adjustment
        wlen_ratio = (self.wlen_grid / seeing_wlen).si

        # Initialize acceptance ratio grid
        corrections = np.empty(
            (self.spec_file.num_exposures, self.num_offset_targets, self.wlen_grid_steps, self.steps_per_exposure),
            dtype=float)

        guided_centroids = []

        # Loop over exposures
        for exp_index in range(self.spec_file.num_exposures):

            # Estimate the actual offset target paths during the exposure
            guided_x, guided_y = self.get_exp_centroids(exp_index)
            guided_centroids.append((guided_x, guided_y))

            # Calculate centroid offsets for each offset target, relative to its nominal fiber center.
            offset = np.sqrt(
                (guided_x - self.offset_x0[:, np.newaxis, np.newaxis])**2 +
                (guided_y - self.offset_y0[:, np.newaxis, np.newaxis])**2)
            
            # Calculate centroid offsets for each offset target, relative to where its fiber center would
            # be if it were designed for the same wavelength as the standard stars.
            offset_std = np.sqrt(
                (guided_x - self.offset_x0_std[:, np.newaxis, np.newaxis])**2 +
                (guided_y - self.offset_y0_std[:, np.newaxis, np.newaxis])**2)

            seeing = self.seeing[exp_index]
            # psf = sdss_25m.get_atmospheric_psf(seeing_wlen, seeing, gauss=False)
            # acceptance_model = sdss_25m.calculate_fiber_acceptance(psf)

            seeing_wlen_adjusted = seeing*wlen_ratio**(-0.2)
            # acceptance_model_grid = map(tpcorr.acceptance_model.AcceptanceModel, seeing_wlen_adjusted)

            for wlen_index in range(self.wlen_grid_steps):
                # Build acceptance model for this wavelength
                acceptance_model = tpcorr.acceptance_model.AcceptanceModel(seeing_wlen_adjusted[wlen_index])

                # Calculate the acceptance fractions for both sets of centroid offsets.
                acceptance = acceptance_model((offset[:,wlen_index,:] / self.pointing.platescale).to(u.arcsec))
                acceptance_std = acceptance_model((offset_std[:,wlen_index,:] / self.pointing.platescale).to(u.arcsec))
                
                # Calculate the acceptance fraction ratios, tabulated for each offset target, wavelength and time.
                # The ratio calculated this way gives the correction of eqn (13).
                corrections[exp_index,:,wlen_index,:] = acceptance_std / acceptance

        # Average the correction over each exposure time slice.
        avg_corrections = np.mean(np.mean(corrections, axis=-1), axis=0)

        return corrections, avg_corrections, guided_centroids

if __name__ == '__main__':
    pass