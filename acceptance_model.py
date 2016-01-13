## Fiber Acceptance Model

import numpy as np

import scipy.interpolate

import astropy.units as u

import galsim

class Telescope(object):
    """
    Represents a telescope.
    """
    def __init__(self,diameter=3.80*u.m, obscuration_area_fraction=0.25, 
                 throughput=0.95*0.77, plate_scale=67.40*u.um/u.arcsec):
        self.diameter = diameter
        self.obscuration_area_fraction = obscuration_area_fraction
        self.throughput = throughput
        self.plate_scale = plate_scale
        self.effective_area = np.pi*diameter**2/4.*(1-obscuration_area_fraction)
    def get_optical_psf(self, wavelength):
        #Convert dimensionless lam/D to arcsec units.
        lam_over_diam_arcsec = ((wavelength/self.diameter)*u.rad).to(u.arcsec)
        # Airy requires floats as inputs, not numpy scalars.
        return galsim.Airy(lam_over_diam=float(lam_over_diam_arcsec.value),
            obscuration=float(np.sqrt(self.obscuration_area_fraction)))
    def get_atmospheric_psf(self, wavelength, fwhm5400, gauss=False):
        wlen_ratio = (wavelength/(5400*u.Angstrom)).si
        assert wlen_ratio == wlen_ratio.value,'wavelength has invalid units.'
        fwhm = fwhm5400.to(u.arcsec).value*wlen_ratio**(-0.2)
        # Kolmogorov requires floats as inputs, not numpy scalars.
        if gauss:
            return galsim.Gaussian(fwhm=float(fwhm))
        else:
            return galsim.Kolmogorov(fwhm=float(fwhm))
    def get_psf(self, wavelength, fwhm5400, rms_jitter=0.1*u.arcsec, gauss=False):
        components = [ self.get_atmospheric_psf(wavelength, fwhm5400, gauss=gauss),self.get_optical_psf(wavelength) ]
        # Include a Gaussian pointing jitter, if requested.
        if rms_jitter is not None:
            components.append(galsim.Gaussian(sigma=rms_jitter.to(u.arcsec).value))
        return galsim.Convolve(components)

def calculate_fiber_acceptance(fiber_diameter, psf, sampling=100, max_offset=2):
    """
    Calculate the fiber acceptance fraction versus offset for a specified PSF.
    
    Args:
        fiber_diameter: Diameter of the fiber to use with explicit angular units.
        psf: PSF model to use, assumed to be specified in arcseconds.
        sampling: Sampling to use for the calculation. Higher samplings take longer
            but give more accurate results.
        max_offset: Maximum centroid offset to calculate, as a ratio of the
            fiber diameter.
    
    Returns:
        tuple: Tuple (offset,acceptance) where offset is given as a ratio of fiber
            diameter and acceptance is a fraction from 0-1.
    """
    # Render the PSF to an image with size fiber_diameter by (max_offset+1)*fiber_diameter.
    diam_arcsec = (fiber_diameter.to(u.arcsec)).value
    width = 2*sampling+1
    height = int(np.ceil((max_offset+1)*width))
    image = galsim.Image(width,height,scale=diam_arcsec/width)
    psf.shift(dx=0.,dy=-0.5*diam_arcsec*max_offset).drawImage(image=image)
    # Prepare a boolean mask of pixels falling inside the fiber aperture.
    xy = np.arange(width) - 0.5*(width-1)
    x,y = np.meshgrid(xy,xy)
    mask = (x**2 + y**2 < (0.5*width)**2)
    # Loop over centroid offsets.
    offset = np.arange(height-width+1)/float(width)
    acceptance = np.empty_like(offset)
    for dy in range(height-width):
        acceptance[dy] = np.sum(image.array[dy:dy+width]*mask)
    #return offset,acceptance
    return scipy.interpolate.interp1d(offset, acceptance, kind='linear', copy=True, bounds_error=True)

class AcceptanceModel(object):
    def __init__(self, seeing_fwhm, fiber_diameter=2*u.arcsec, sampling=100, max_offset=0.75):
        """
        Calculate the fiber acceptance fraction versus offset for a specified PSF.

        Args:
            fiber_diameter: Diameter of the fiber to use with explicit angular units.
            psf: PSF model to use, assumed to be specified in arcseconds.
            sampling: Sampling to use for the calculation. Higher samplings take longer
                but give more accurate results.
            max_offset: Maximum centroid offset to calculate, as a ratio of the
                fiber diameter.
        """
        # Render the PSF to an image with size fiber_diameter by (max_offset+1)*fiber_diameter.
        diam_arcsec = (fiber_diameter.to(u.arcsec)).value
        width = 2 * sampling + 1
        pixel_scale = diam_arcsec / width
        height = int(np.ceil((max_offset + 1) * width))
        sigma_arcsec = seeing_fwhm.to(u.arcsec).value / (2. * np.sqrt(2. * np.log(2)))
        xx = pixel_scale * (np.arange(width) - 0.5 * (width - 1))
        yy = pixel_scale * (np.arange(height) - 0.5 * (width - 1))
        x, y = np.meshgrid(xx, yy, sparse=True, copy=False)
        image = np.exp(-(x**2 + y**2) / (2. * sigma_arcsec**2)) * pixel_scale**2 / (2 * np.pi * sigma_arcsec**2)
        # Prepare a boolean mask of pixels falling inside the fiber aperture.
        mask = xx**2 + xx[:, np.newaxis]**2 < (0.5 * diam_arcsec)**2
        # Loop over centroid offsets.
        offset = np.arange(height - width + 1) * pixel_scale
        acceptance = np.zeros_like(offset)
        for dy in range(height - width + 1):
            acceptance[dy] = np.sum(image[dy:dy + width] * mask)
        ##plt.plot(offset, acceptance)
        self.interpolation = scipy.interpolate.interp1d(
            offset, acceptance, kind='linear', copy=True, bounds_error=True)

    def __call__(self, centroid_offsets):
        return self.interpolation(centroid_offsets.to(u.arcsec).value)

if __name__ == '__main__':
    pass
