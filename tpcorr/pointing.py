## Atmospheric Refraction Model

import numpy as np

import scipy.interpolate

import astropy.time
import astropy.coordinates
import astropy.units as u

import specsim

import tpcorr.distortion

# helper function for shifting angles to desired range
def normalize_angle(angle):
    while angle <= -180:
        angle += 360
    while angle > 180:
        angle -= 360
    return angle

class Pointing(object):
    """Represents an observer pointing at a fixed sky location.
    """
    def __init__(self, ra_center, dec_center):
        self.plate_center = astropy.coordinates.ICRS(ra=ra_center, dec=dec_center)
        self.plate_fid = astropy.coordinates.ICRS(ra=ra_center, dec=dec_center + 1.5 * u.deg)
        self.platescale = 217.7358 * u.mm / u.deg
        self.wlen0 = 5500 * u.Angstrom
        self.relhum = 0.2
        try:
            self.where = astropy.coordinates.EarthLocation.of_site('apo')
        except astropy.coordinates.errors.UnknownSiteException:
            self.where = astropy.coordinates.EarthLocation(lat=32.7797556*u.deg, lon=-(105+49./60.+13/3600.)*u.deg, height=2797*u.m)
        self.distortion_model = tpcorr.distortion.get_optical_distortion_model(330.0 * u.mm, self.platescale)
        self.chromatic_model = tpcorr.distortion.get_chromatic_distortion_model(self.platescale)

    def transform(self, targets, tai, wlen, temperature, pressure, do_chromatic_distortion=True, extrap_wlen=False):
        """Transform from sky coordinates to focal plane coordinates.
        
        Args:
            targets: astropy.coordinates.SkyCoord
            tai: float or numpy.ndarray
            wlen: astropy.units.quantity.Quantity
            temperature: astropy.units.quantity.Quantity
            pressure: astropy.units.quantity.Quantity or None to calculate based
                on temperature and elevation
            
        Returns:
            tuple x, y of astropy.units.quantity.Quantity objects giving focal plane
            positions in length units, broadcast over the input args.
        """

        # Initialize time objects from the input TAI values in MJD seconds.
        when = astropy.time.Time(tai/86400., format='mjd', scale='tai', location=self.where)
                
        # Calculate the Alt-Az path of the telescope boresight at the plate design wavelength (5500A).
        obs_model0 = specsim.transform.create_observing_model(
                    where=self.where, when=when, wavelength=self.wlen0, pressure=pressure,
                    temperature=temperature, relative_humidity=self.relhum)
        altaz0 = specsim.transform.sky_to_altaz(self.plate_center, obs_model0)
        alt0, az0 = altaz0.alt, altaz0.az
    
        # Calculate the Alt-Az paths of each target over the input wavelength grid.
        obs_model = specsim.transform.create_observing_model(
            where=self.where, when=when, wavelength=wlen, pressure=pressure,
            temperature=temperature, relative_humidity=self.relhum)
        altaz = specsim.transform.sky_to_altaz(targets, obs_model)

        # Convert each target's Alt-Az into local X, Y focal plane coordinates.
        x, y = specsim.transform.altaz_to_focalplane(
            altaz.alt, altaz.az, altaz0.alt, altaz0.az, self.platescale)
        # Flip y to match the handedness of the XFOCAL, YFOCAL coordinate system.
        y = -y

        # Rotate the focal plane so that +y points towards a point that is offset from
        # the plate center along DEC by +1.5 degrees.
        altaz_fid = specsim.transform.sky_to_altaz(self.plate_fid, obs_model0)
        x_fid, y_fid = specsim.transform.altaz_to_focalplane(
            altaz_fid.alt, altaz_fid.az, altaz0.alt, altaz0.az, self.platescale)
        angle = np.arctan2(x_fid.si, -y_fid.si)
        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)
        x_rot = x * cos_angle - y * sin_angle
        y_rot = x * sin_angle + y * cos_angle

        # Apply radial optical distortions.
        r = np.sqrt(x_rot**2 + y_rot**2)
        distortion = ((r + self.distortion_model(r)) / r).si

        x_dist = distortion * x_rot
        y_dist = distortion * y_rot

        if do_chromatic_distortion:
            r_dist = np.sqrt(x_dist**2 + y_dist**2)
            dr5000 = self.chromatic_model(r_dist, 5000, extrap_wlen)
            dr = np.empty_like(r_dist)
            # ugh broadcasting...
            if wlen.isscalar:
                dr = self.chromatic_model(r_dist, wlen, extrap_wlen) - dr5000
            elif len(wlen.shape) == 2:
                for iw,w in enumerate(wlen.flatten()):
                    dr[:,iw] = self.chromatic_model(r_dist[:,iw], w.to(u.Angstrom).value, extrap_wlen) - dr5000[:,iw]
            chromatic_distortion = ((r_dist + dr) / r_dist).si

            x_dist = chromatic_distortion * x_dist
            y_dist = chromatic_distortion * y_dist
        
        return x_dist, y_dist, altaz.alt, altaz.az
    
    def hour_angle(self, tai):
        """Convert TAI to the hour angle of this plate's RA"""
        when = astropy.time.Time(tai/86400., format='mjd', scale='tai', location=self.where)
        return when.sidereal_time('apparent') - self.plate_center.ra

if __name__ == '__main__':
    p = Pointing(ra_center=180*u.deg, dec_center=25*u.deg)
