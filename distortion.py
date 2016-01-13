## Optical Distortion Model

# The chromatic distortion is described [here](https://trac.sdss3.org/tracmailman/browser/private/sdss3-infrastructure/all/505.html) 
# with corrected tables appearing [here](https://trac.sdss3.org/tracmailman/browser/private/sdss3-infrastructure/all/946.html). 
# The table of *principal ray heights* are tabulated in [platedesign svn repo](https://trac.sdss.org/browser/repo/sdss/platedesign/trunk/data/sdss/image-heights.txt) 
# and used in [this IDL code](https://trac.sdss3.org/browser/repo/platedesign/trunk/pro/plate/apo_rdistort.pro). 
# Note that the relative distortion number only have one significant digit, so are of limited use.

import numpy as np
import matplotlib.pyplot as plt

import scipy.interpolate

import astropy.units as u

# https://trac.sdss.org/browser/repo/sdss/platedesign/trunk/data/sdss/plParam.par
optical_distortion_coefs = np.array([
        -0.000137627, -0.00125238, 1.5447e-09, 8.23673e-08, -2.74584e-13, -1.53239e-12, 6.04194e-18,
        1.38033e-17, -2.97064e-23, -3.58767e-23, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

def get_optical_distortion_model(r_max, platescale, npts=100):
    r = np.linspace(0, r_max.to(u.mm).value, npts)*u.mm
    asin_r = np.arcsin((r / platescale).to(u.rad).value) * u.rad
    dr = np.polynomial.polynomial.polyval((asin_r * platescale).to(u.mm).value, optical_distortion_coefs)
    interpolator = scipy.interpolate.interp1d(asin_r.value, dr, copy=True, kind='cubic')
    def model(r):
        # Input r should be a Quantity with length units.
        return interpolator(np.arcsin((r / platescale).to(u.rad).value)) * u.mm
    return model

def plot_optical_distortion(r_max=325*u.mm, platescale=217.7358*u.mm/u.deg):
    model = get_optical_distortion_model(r_max, platescale)
    r = np.linspace(0, r_max.to(u.mm).value, 100)*u.mm
    dr = (model(r) / platescale).to(u.arcsec)
    plt.plot(r, dr)
    plt.grid()
    plt.xlabel('Undistorted radius [mm]')
    plt.ylabel('Radial distortion at 5000A [arcsec]')


distortion_wlen = np.array([
    4000., 5300., 5500., 6000., 8000., 10000., 15350., 15950., 16550. ]) * u.Angstrom

distortions_at_5300 = np.array([
    -0.00,  36.26,  72.53, 108.84, 145.18, 181.53, 217.90, 254.29, 290.77, 327.44
]) * u.mm

distortions_relative_to_5300 = np.array([
    [  0.000, -0.002, -0.003, -0.004, -0.005, -0.005, -0.005, -0.004, -0.002,  0.003 ],
    # Added row of zeros for 5300A
    [  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000 ],
    [ -0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.000 ],
    [  0.000,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001, -0.001 ],
    [  0.000,  0.001,  0.003,  0.003,  0.004,  0.004,  0.004,  0.003,  0.002, -0.003 ],
    [  0.000,  0.002,  0.004,  0.005,  0.005,  0.005,  0.005,  0.005,  0.003, -0.004 ],
    [ -0.000,  0.003,  0.006,  0.007,  0.008,  0.008,  0.008,  0.008,  0.004, -0.006 ],
    [ -0.000,  0.003,  0.006,  0.008,  0.008,  0.009,  0.009,  0.008,  0.004, -0.006 ],
    [  0.000,  0.004,  0.006,  0.008,  0.009,  0.009,  0.009,  0.008,  0.004, -0.007 ]
]) * u.mm


def get_chromatic_distortion_model(platescale):
    """
    https://trac.sdss3.org/browser/repo/platedesign/trunk/pro/plate/apo_rdistort.pro
    """
    # Chromatic distortions are tabulated at 10 radii, corresponding to
    # sin(r) = 0', 10', ..., 90' (with arcmins converted to radians).
    sin_r_table = (np.arange(10) * 10 * u.arcmin).to(u.rad)
    # Calculate the corresponding radii in mm.
    r_table = (np.arcsin(sin_r_table.value) * u.rad * platescale).to(u.mm)
    
    # Calculate fractional distortions relative to 5300A
    d5300 = distortions_relative_to_5300 / distortions_at_5300 
    # We are copying the original IDL code here, but setting these to zero might
    # make more sense.
    d5300[:, 0] = d5300[:, 1]
    
    # Calculate additive distortions relative to 5500A.
    assert distortion_wlen[2] == 5500. * u.Angstrom
    d5500 = d5300 - d5300[2]
    
    # Build a linear interpolator in wavelength of additive distortions relative to 5500A.
    wlen_interpolator = scipy.interpolate.interp1d(distortion_wlen, d5500, axis=0,
                                                   kind='linear', copy=True, bounds_error=True)
    
    def model(r, wlen):
        # Clip wavelengths to the tabulated limits.
        wlen = np.clip(wlen, distortion_wlen[0], distortion_wlen[-1])
        # Calculate the additive relative distortion at each wavelength.
        radial_distortion = wlen_interpolator(wlen)
        r_interpolator = scipy.interpolate.interp1d(r_table, radial_distortion, kind='cubic',
                                                    copy=False, bounds_error=True)
        return r * r_interpolator(r)
    
    return model


def plot_chromatic_distortion(wlen=4000*u.Angstrom, r_max=325*u.mm, platescale=217.7358*u.mm/u.deg):
    model = get_chromatic_distortion_model(platescale)
    r = np.linspace(0, r_max.to(u.mm).value, 100)*u.mm
    dr = (model(r, wlen.to(u.Angstrom).value) / platescale).to(u.arcsec)
    plt.plot(r, dr)
    plt.grid(True)
    plt.xlabel('Undistorted radius [mm]')
    plt.ylabel('Radial distortion relative to 5500A [arcsec]')


if __name__ == '__main__':
    plot_optical_distortion()
    plt.show()
    wlen_array = np.linspace(3500, 10500, 15)*u.Angstrom
    for wlen in wlen_array:
        plot_chromatic_distortion(wlen=wlen)
    plt.show()
