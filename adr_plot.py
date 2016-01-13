#!/usr/bin/env python

# plots differntial refraction relative to 4000A and 5400A (adr_plot.pdf)

import numpy as np
import matplotlib as mpl
mpl.rcParams.update({'font.size': 18})
mpl.rcParams.update({'savefig.dpi': 200})
mpl.rcParams.update({'savefig.bbox': 'tight'})

import matplotlib.pyplot as plt

from astropy.utils.data import download_file
from astropy.utils import iers
iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL, cache=True))

import astropy.units as u
import astropy.coordinates
import astropy.time


def main():
    # observatory 
    apo = astropy.coordinates.EarthLocation.of_site('apo')
    pressure = 72.80555*u.kPa
    temperature = 15*u.deg_C
    when = astropy.time.Time.now()

    # pointings
    num_alt_steps = 5
    alt = np.linspace(50, 90, num_alt_steps)*u.deg
    az = np.zeros(num_alt_steps)*u.deg
    pointings = astropy.coordinates.SkyCoord(az, alt, frame='altaz', location=apo, obstime=when)

    # wavelength array
    num_wlen_steps = 71
    wlen_array = np.linspace(3500, 10500, num_wlen_steps)*u.Angstrom

    # calculate refracted altitudes at each alt, wavelength
    refracted_alt_grid = np.empty((num_wlen_steps, num_alt_steps))*u.deg
    for i,wlen in enumerate(wlen_array):
        refracted_alt_grid[i] = pointings.transform_to(
            astropy.coordinates.AltAz(location=apo, obswl=wlen, pressure=pressure, temperature=temperature)).alt 

    # create figure
    fig = plt.figure(figsize=(8,6))

    # plot differential refraction (relative to 4000A and 5400A) at each altitude
    for i in range(num_alt_steps):
        i4000, i5400 = 5, 19
        # change line alpha with altitude
        alpha = 1 - (i * 0.2)
        plt.plot(wlen_array, (refracted_alt_grid[:,i]-refracted_alt_grid[i4000,i]).to(u.arcsec), c='b', alpha=alpha)
        plt.plot(wlen_array, (refracted_alt_grid[i5400,i]-refracted_alt_grid[:,i]).to(u.arcsec), c='r', alpha=alpha)

        # add altitude labels
        i3600, i5800 = 1, 23
        plt.text(wlen_array[i3600].value, ((refracted_alt_grid[i5400,i]-refracted_alt_grid[i3600,i]).to(u.arcsec)-0.075*u.arcsec).value, r'$%d^\circ$'%alt[i].value, fontsize=14, color='r')
        plt.text(wlen_array[i5800].value, ((refracted_alt_grid[i5800,i]-refracted_alt_grid[i4000,i]).to(u.arcsec)-0.085*u.arcsec).value, r'$%d^\circ$'%alt[i].value, fontsize=14, color='b')
    # plot labels and limits
    plt.ylabel('ADR $(\mathrm{arcseconds})$')
    plt.xlabel('Observed Wavelength $\lambda$ $(\AA)$')
    plt.xlim(3500, 6000)
    plt.ylim(-1, 0.5)

    # vertical dashed lines
    crossover_wlen = 4539*u.Angstrom
    plt.axvline(4000, ls='--', c='b', zorder=0)
    plt.axvline(5400, ls='--', c='r', zorder=0)
    plt.axvline(crossover_wlen.value, ls='-', c='k', zorder=0)

    # horizontal grid lines
    plt.grid(axis='y')

    # save figure
    plt.savefig('adr_plot.pdf')

if __name__ == '__main__':
    main()