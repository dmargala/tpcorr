#!/usr/bin/env python

import argparse

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.size': 18})
mpl.rcParams.update({'savefig.dpi': 200})
mpl.rcParams.update({'savefig.bbox': 'tight'})

import matplotlib.pyplot as plt

import scipy.interpolate

import astropy.units as u

import tpcorr

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--wlen', type=float, default=5400.,
        help='Observing wavelength')
    args = parser.parse_args()

    # Observing wavelength
    wlen = args.wlen * u.Angstrom

    # Initialize telescope model
    sdss_25m = tpcorr.acceptance_model.Telescope()

    # psf fwhm 
    fwhm0 = 1.5 * u.arcsec
    fwhm_array = np.array([1.2, 1.5, 1.8]) * u.arcsec
    linestyles = ['-.', 'solid', 'dashed']

    # offset sampling grid
    dmin, dmax, nd = (0, 1.5, 50)
    offsets = np.linspace(dmin, dmax, nd)
    offset_std_grid, offset_grid = np.meshgrid(offsets, offsets)
    edges = np.linspace(dmin, dmax, nd + 1)

    # Levels for acceptance ratio contours
    acceptance_ratio_levels = np.linspace(0.5, 0.9, 3)
    # acceptance_ratio_levels_inverse = [1.0/l for l in acceptance_ratio_levels[::-1]]
    
    plt.figure(figsize=(8,6))
    for i,(fwhm,linestyle) in enumerate(zip(fwhm_array,linestyles)):
        # Draw reference contours for Gaussian atmosphere
        psf = sdss_25m.get_atmospheric_psf(wlen, fwhm, gauss=True)
        acceptance = sdss_25m.calculate_fiber_acceptance(psf)
        # acceptance = tpcorr.acceptance_model.AcceptanceModel(fwhm)
        acceptance_ratio = acceptance(offset_std_grid) / acceptance(offset_grid)
        contours_ratio = plt.contour(offset_grid, offset_std_grid, acceptance_ratio, acceptance_ratio_levels,
            colors='green', linewidths=1, linestyles=linestyle)
        contours_ratio = plt.contour(offset_grid, offset_std_grid, acceptance_ratio, 1.0/acceptance_ratio_levels,
            colors='green', linewidths=1, linestyles=linestyle)

    # Draw reference contours for Kolmogorov atmosphere, use the reference fwhm
    psf = sdss_25m.get_atmospheric_psf(wlen, fwhm0, gauss=False)
    acceptance = sdss_25m.calculate_fiber_acceptance(psf)
    acceptance_ratio = acceptance(offset_std_grid) / acceptance(offset_grid)
    contours_ratio = plt.contour(offset_grid, offset_std_grid, acceptance_ratio, acceptance_ratio_levels,
        colors='black', linewidths=1, linestyles='solid')
    # Add contour labels
    plt.clabel(contours_ratio, fontsize=11, fmt=lambda l: '%.2f'%l)
    contours_ratio = plt.contour(offset_grid, offset_std_grid, acceptance_ratio, 1.0/acceptance_ratio_levels,
        colors='black', linewidths=1, linestyles='solid')
    # Add contour labels
    plt.clabel(contours_ratio, fontsize=11, fmt=lambda l: '%.2f'%l)

    # draw 1 to 1 line 
    plt.plot(offsets, offsets, color='green', ls='-', lw=1)
        
    # Set aspect ratio and plot limits
    plt.gca().set_aspect('equal')
    plt.xlim(dmin, dmax)
    plt.ylim(dmin, dmax)

    # Add title and axis labels
    plt.xlabel(r'$d_i^\ast(\lambda,\lambda_i,h_\mathrm{obs})$ $(\mathrm{arcseconds})$')
    plt.ylabel(r'$d_i^\ast(\lambda,\lambda_c,h_\mathrm{obs})$ $(\mathrm{arcseconds})$')
    plt.title(r'$\lambda$ = %d $\AA$' % wlen.value)

    # Save figure
    plt.savefig('acceptance_contours.pdf')


if __name__ == '__main__':
    main()
