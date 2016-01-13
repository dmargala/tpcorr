#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
mpl.rcParams.update({'font.size': 18})
mpl.rcParams.update({'savefig.dpi': 200})
mpl.rcParams.update({'savefig.bbox': 'tight'})

import matplotlib.pyplot as plt

import scipy.interpolate

import astropy.units as u

from acceptance_model import *



def plot_offset_acceptance(t, D=2*u.arcsec, wlen=5400*u.Angstrom, fwhm1=1.2, fwhm2=1.8, nfwhm=3, 
                           noffset=50, offset_min=0, offset_max=2):
    cmap = mpl.cm.ScalarMappable(
        mpl.colors.Normalize(vmin=fwhm1-.1,vmax=fwhm2+.1),'spectral')
    
    plt.subplot(1,1,1)
    fwhm_vec = np.linspace(fwhm1,fwhm2,nfwhm)
    A_vec = np.empty((nfwhm, noffset,))
    A_gauss_vec = np.empty((nfwhm, noffset,))
    A_atmos_vec = np.empty((nfwhm, noffset,))
    offsets_vec = np.linspace(offset_min, offset_max,noffset)
    
    for i,fwhm in enumerate(fwhm_vec):
        psf = t.get_psf(wlen,fwhm*u.arcsec)
        interpolator = calculate_fiber_acceptance(D,psf)
        for j,offset in enumerate(offsets_vec):
            A_vec[i,j] = interpolator((offset*u.arcsec/D).si.value)
            
        psf_atmos = t.get_atmospheric_psf(wlen,fwhm*u.arcsec)
        interpolator = calculate_fiber_acceptance(D, psf_atmos)
        for j,offset in enumerate(offsets_vec):
            A_atmos_vec[i,j] = interpolator((offset*u.arcsec/D).si.value)
                    
        psf_atmos_gauss = t.get_atmospheric_psf(wlen,fwhm*u.arcsec, gauss=True)
        interpolator = calculate_fiber_acceptance(D, psf_atmos_gauss)
        for j,offset in enumerate(offsets_vec):
            A_gauss_vec[i,j] = interpolator((offset*u.arcsec/D).si.value)
                    
    linestyles = ['-.','-','--']
    for i,fwhm in enumerate(fwhm_vec):
        color = cmap.to_rgba(fwhm)
        # plt.plot(offsets_vec,A_vec[i],color=color, ls='-.',)
        plt.plot(offsets_vec,A_atmos_vec[i], color='black', lw=2, ls=linestyles[i], label=(r'%.1f${}^{\prime\prime}$'%fwhm))
        plt.plot(offsets_vec,A_gauss_vec[i], color='green', lw=2, ls=linestyles[i])
        
    plt.xlabel(r'Centroid Offset $d$ $(\mathrm{arcseconds})$')
    plt.ylabel(r'Fiber Acceptance $A(d)$')

def main():

    sdss_25m = Telescope(diameter=2.5*u.m, obscuration_area_fraction=0.27, plate_scale=217.7358*u.mm/u.deg)

    plt.figure(figsize=(8,6))
    plot_offset_acceptance(sdss_25m, fwhm1=1.2, fwhm2=1.8, offset_max=2, noffset=200)
    plt.axvline(1, color='gray', ls='-')
    plt.xlim(0, 1.5)
    plt.ylim(0, 1)
    plt.grid(True)
    plt.savefig('acceptance_plot.pdf')

if __name__ == '__main__':
    main()