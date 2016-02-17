## Fiber Acceptance Model

import numpy as np

import scipy.interpolate

import astropy.units as u

class AcceptanceModel(object):
    def __init__(self, seeing_fwhm, fiber_diameter=2*u.arcsec, sampling=100, max_offset=2.0):
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

class Telescope(object):
    """
    Represents a telescope.
    """
    def __init__(self, diameter=2.5*u.m, obscuration_area_fraction=0.27, plate_scale=217.7358*u.mm/u.deg):
        import galsim
        self.galsim = galsim
        self.diameter = diameter
        self.obscuration_area_fraction = obscuration_area_fraction
        self.plate_scale = plate_scale
        self.effective_area = np.pi*diameter**2/4.*(1-obscuration_area_fraction)
    def get_optical_psf(self, wavelength):
        #Convert dimensionless lam/D to arcsec units.
        lam_over_diam_arcsec = ((wavelength/self.diameter)*u.rad).to(u.arcsec)
        # Airy requires floats as inputs, not numpy scalars.
        return self.galsim.Airy(lam_over_diam=float(lam_over_diam_arcsec.value),
            obscuration=float(np.sqrt(self.obscuration_area_fraction)))
    def get_atmospheric_psf(self, wavelength, fwhm5400, gauss=False):
        wlen_ratio = (wavelength/(5400*u.Angstrom)).si
        assert wlen_ratio == wlen_ratio.value,'wavelength has invalid units.'
        fwhm = fwhm5400.to(u.arcsec).value*wlen_ratio**(-0.2)
        # Kolmogorov requires floats as inputs, not numpy scalars.
        if gauss:
            return self.galsim.Gaussian(fwhm=float(fwhm))
        else:
            return self.galsim.Kolmogorov(fwhm=float(fwhm))
    def get_psf(self, wavelength, fwhm5400, rms_jitter=0.1*u.arcsec, gauss=False):
        components = [ self.get_atmospheric_psf(wavelength, fwhm5400, gauss=gauss),self.get_optical_psf(wavelength) ]
        # Include a Gaussian pointing jitter, if requested.
        if rms_jitter is not None:
            components.append(self.galsim.Gaussian(sigma=rms_jitter.to(u.arcsec).value))
        return self.galsim.Convolve(components)

    def calculate_fiber_acceptance(self, psf, fiber_diameter=2*u.arcsec, sampling=100, max_offset=2, return_arrays=False):
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
        image = self.galsim.Image(width,height,scale=diam_arcsec/width)
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
        if return_arrays:
            return offset*fiber_diameter.to(u.arcsec).value, acceptance
        else:
            return scipy.interpolate.interp1d(
                offset*fiber_diameter.to(u.arcsec).value, acceptance, kind='linear', copy=True, bounds_error=True)

if __name__ == '__main__':

    std_fwhm = 1.5 * u.arcsec
    std_wlen = 5400 * u.Angstrom

    acceptance_model = AcceptanceModel(std_fwhm)

    wlen_array = np.linspace(3500., 10500., 15) * u.Angstrom
    wlen_ratio = (wlen_array / std_wlen).si

    offset_grid = np.zeros((5, 15))
    offset_grid += np.linspace(0, 1., 15)[np.newaxis, :]
    offset_grid += np.linspace(0,.05, 5)[:, np.newaxis]

    acceptance_grid = np.empty_like(offset_grid)

    fwhm_array = std_fwhm*wlen_ratio**(-0.2)
    amodel_array = map(AcceptanceModel, fwhm_array)

    for i, wlen in enumerate(wlen_array):
        amodel = amodel_array[i]
        acceptance_grid[:,i] = amodel(offset_grid[:,i]*u.arcsec)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8,6))
    plt.plot(wlen_array, fwhm_array)
    plt.axvline(std_wlen.value, c='k', ls='--')
    plt.axhline(std_fwhm.value, c='k', ls='--')
    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('PSF FWHM (arcsec)')
    plt.xlim(wlen_array[0].value, wlen_array[-1].value)
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(8,6))
    for i in range(5):
        plt.plot(wlen_array, acceptance_grid[i])
    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('Acceptance')
    plt.xlim(wlen_array[0].value, wlen_array[-1].value)
    plt.grid(True)
    
    print 'Saving figure: acceptance_model_test.png'
    plt.savefig('acceptance_model_test.png')





