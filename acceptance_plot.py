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
        mpl.colors.Normalize(vmin=fwhm1-.1, vmax=fwhm2+.1),'spectral')
    
    plt.subplot(1,1,1)
    fwhm_vec = np.linspace(fwhm1, fwhm2, nfwhm)
    A_vec = np.empty((nfwhm, noffset,))
    A_gauss_vec = np.empty((nfwhm, noffset,))
    A_atmos_vec = np.empty((nfwhm, noffset,))
    offsets_vec = np.linspace(offset_min, offset_max, noffset)
    
    for i,fwhm in enumerate(fwhm_vec):
        psf = t.get_psf(wlen, fwhm*u.arcsec)
        interpolator = calculate_fiber_acceptance(D, psf)
        for j,offset in enumerate(offsets_vec):
            A_vec[i,j] = interpolator((offset*u.arcsec/D).si.value)
            
        psf_atmos = t.get_atmospheric_psf(wlen, fwhm*u.arcsec)
        interpolator = calculate_fiber_acceptance(D, psf_atmos)
        for j,offset in enumerate(offsets_vec):
            A_atmos_vec[i,j] = interpolator((offset*u.arcsec/D).si.value)
                    
        psf_atmos_gauss = t.get_atmospheric_psf(wlen, fwhm*u.arcsec, gauss=True)
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

def plot_offset_acceptance_ratio(t, D=2*u.arcsec, wlen=5400*u.Angstrom, fwhm=1.5, sampling=100):

    plt.subplot(1,1,1)

    psf = t.get_atmospheric_psf(wlen,fwhm*u.arcsec)
    offsets, A = calculate_fiber_acceptance(D, psf, return_arrays=True)
    A_interp = scipy.interpolate.interp1d(offsets, A)
            
    gauss_psf = t.get_atmospheric_psf(wlen, fwhm*u.arcsec,gauss=True)
    offsets, A_gauss = calculate_fiber_acceptance(D, gauss_psf, return_arrays=True)
    A_gauss_interp = scipy.interpolate.interp1d(offsets, A_gauss)
    
    offset_vec = np.arange(.1, .6, .1)
    cmap = mpl.cm.ScalarMappable(
        mpl.colors.Normalize(vmin=-.1, vmax=.6),'spectral')
    for i,d in enumerate(offset_vec):   
        color = cmap.to_rgba(d)
        plt.plot(offsets*D.to(u.arcsec).value, A/A_interp(d/D.to(u.arcsec).value), color=color, label=r'%.1f${}^{\prime\prime}$'%d)
        plt.plot(offsets*D.to(u.arcsec).value, A_gauss/A_gauss_interp(d/D.to(u.arcsec).value), color=color, ls='--')

    plt.xlim(0, 2)
    plt.xlabel(r'Centroid Offset $d_{4000}$ $[\ {}^{\prime\prime}]\ $')
    plt.ylabel(r'$A(\sigma_\mathrm{PSF},\ d_{4000}) / A(\sigma_\mathrm{PSF},\ d_{5400})$')

def plot_offset_acceptance_ratio_ratio(t, D=2*u.arcsec, wlen=5400*u.Angstrom, fwhm=1.5, sampling=100):

    plt.subplot(1,1,1)

    psf = t.get_atmospheric_psf(wlen, fwhm*u.arcsec)
    offsets, A = calculate_fiber_acceptance(D, psf, return_arrays=True)
    A_interp = scipy.interpolate.interp1d(offsets, A)
            
    gauss_psf = t.get_atmospheric_psf(wlen, fwhm*u.arcsec, gauss=True)
    offsets, A_gauss = calculate_fiber_acceptance(D, gauss_psf, return_arrays=True)
    A_gauss_interp = scipy.interpolate.interp1d(offsets, A_gauss)
    
    offset_vec = np.arange(0, 1.3, 0.3)
    cmap = mpl.cm.ScalarMappable(
        mpl.colors.Normalize(vmin=-.2, vmax=1.4),'spectral')
    for i,d in enumerate(offset_vec):   
        color = cmap.to_rgba(d)
        plt.plot(offsets*D.to(u.arcsec).value, (A_gauss/A_gauss_interp(d/D.to(u.arcsec).value))/(A/A_interp(d/D.to(u.arcsec).value)), color=color, label=r'%.1f${}^{\prime\prime}$'%d)

    plt.xlim(0, 2)
    plt.xlabel(r'Centroid Offset $d_{4000}$ $[\ {}^{\prime\prime}]\ $')
    plt.ylabel(r'$C_{Gaussian} / C_{Kolmogorov}$')

def plot_wlen_acceptance(t, wlen_vec=[3500,5400,10000], D=2*u.arcsec, fwhm=1.5,
                           noffset=50, offset_min=0, offset_max=2):
    plt.subplot(1,1,1)
    A_vec = np.empty((len(wlen_vec), noffset,))
    A_gauss_vec = np.empty((len(wlen_vec), noffset,))
    offsets_vec = np.linspace(offset_min, offset_max, noffset)
    
    for i,wlen in enumerate(wlen_vec):
        psf = t.get_atmospheric_psf(wlen*u.Angstrom, fwhm*u.arcsec)
        offsets, acceptance = calculate_fiber_acceptance(D, psf, return_arrays=True)
        interpolator = scipy.interpolate.interp1d(offsets, acceptance)
        for j,offset in enumerate(offsets_vec):
            A_vec[i,j] = interpolator((offset*u.arcsec/D).si.value)
            
        gauss_psf = t.get_atmospheric_psf(wlen*u.Angstrom, fwhm*u.arcsec, gauss=True)
        offsets, acceptance = calculate_fiber_acceptance(D, gauss_psf, return_arrays=True)
        interpolator = scipy.interpolate.interp1d(offsets, acceptance)
        for j,offset in enumerate(offsets_vec):
            A_gauss_vec[i,j] = interpolator((offset*u.arcsec/D).si.value)
                    

    linestyles = [':','-','--']
    for i,wlen in enumerate(wlen_vec):
        plt.plot(offsets_vec, A_vec[i], lw=2, ls=linestyles[i], label=('%d'%wlen), color='black')
        plt.plot(offsets_vec, A_gauss_vec[i], lw=2, ls=linestyles[i], color='red')
        
    plt.legend(title='Wavelength $(\AA)$', fontsize=14)
    plt.xlabel(r'Centroid Offset $d$ $(\mathrm{arcseconds})$')
    plt.ylabel(r'Fiber Acceptance $A$')
    plt.xlim(0, 2)


def main():

    sdss_25m = Telescope(diameter=2.5*u.m, obscuration_area_fraction=0.27, plate_scale=217.7358*u.mm/u.deg)

    plt.figure(figsize=(8,6))
    plot_offset_acceptance(sdss_25m, fwhm1=1.2, fwhm2=1.8, offset_max=2, noffset=200)
    plt.axvline(1, color='gray', ls='-')
    plt.xlim(0, 1.5)
    plt.ylim(0, 1)
    plt.grid(True)
    plt.savefig('acceptance_plot.pdf')

    plt.figure(figsize=(8,6))
    plot_offset_acceptance_ratio(sdss_25m)
    plt.legend(title='$d_{5400}$', fontsize=14)
    plt.axvline(1, color='gray', ls='-')
    plt.grid(True)
    plt.savefig('acceptance_ratio_plot.pdf')

    plt.figure(figsize=(8,6))
    plot_offset_acceptance_ratio_ratio(sdss_25m)
    plt.axvline(1, color='gray', ls='-')
    plt.legend(title='$d_{5400}$', fontsize=14)
    plt.ylim(.6,1.2)
    plt.grid(True)
    plt.savefig('acceptance_ratio2_plot.pdf')

    plt.figure(figsize=(8,6))
    plot_wlen_acceptance(sdss_25m)
    plt.axvline(1, color='gray', ls='-')
    plt.savefig('acceptance_wlen_plot.pdf')

if __name__ == '__main__':
    main()