#!/usr/bin/env python

import argparse

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.size': 18})
mpl.rcParams.update({'savefig.dpi': 200})
mpl.rcParams.update({'savefig.bbox': 'tight'})

import matplotlib.pyplot as plt

import scipy.stats
import h5py

def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tpcorr', type=str, default=None,
        help='throughput correction filename, required')
    parser.add_argument('--output', type=str, default=None,
        help='output filename')
    args = parser.parse_args()


    tpcorr = h5py.File(args.tpcorr, 'r')

    tpcorr_old = h5py.File('/Users/Daniel/source/qusp/tpcorr-0.1/tpcorr.hdf5','r')
    wave = np.squeeze(tpcorr['wave'].value)

    dha_threshold = 1.25 #1.25

    tpcorrs = []
    opt_tpcorrs = []
    has = []
    alts = []
    design_has = []
    design_alts = []
    plate_mjd_list = []
    for plate in tpcorr.keys():
        if plate == 'wave':
            continue
        for mjd in tpcorr[plate].keys():
            plate_mjd_list.append((plate, mjd))
            grp = tpcorr['{}/{}'.format(plate, mjd)]

            ha, design_ha = grp.attrs['ha'], grp.attrs['design_ha']
            alt, design_alt = grp.attrs['alt'], grp.attrs['design_alt']

            has.append(grp.attrs['ha'])
            design_has.append(grp.attrs['design_ha'])
            alts.append(grp.attrs['alt'])
            design_alts.append(grp.attrs['design_alt'])

            for fiber in grp.keys():
                tpcorrs.append(grp[fiber].value.tolist())
            if abs(design_ha-ha) < dha_threshold:
                for fiber in tpcorr_old[plate][mjd]:
                    opt_tpcorrs.append(tpcorr_old[plate][mjd][fiber].value.tolist())

    num_corrections = len(tpcorrs)
    num_optimal_obs = len(opt_tpcorrs)
    print 'Number of corrections: ', num_corrections

    # Fraction of 
    print float(num_corrections-num_optimal_obs)/num_corrections

    # Convert lists to numpy arrays
    tpcorrs = np.array(tpcorrs)
    opt_tpcorrs = np.array(tpcorrs)
    mean_ha = np.array(has)
    design_ha = np.array(design_has)
    mean_alt = np.array(alts)
    design_alt = np.array(design_alts)

    delta_ha = mean_ha - design_ha
    optimal_obs = np.where(np.abs(delta_ha) < dha_threshold)[0]
    optimal_plate_mjd_list = [plate_mjd for plate_mjd,is_optimal in zip(plate_mjd_list, optimal_obs) if is_optimal]

    # calculate cumulative dist values for a normal distribution
    percentiles = [100*scipy.stats.norm.cdf(sigma) for sigma in [-2,-1,0,1,2]]

    # N-sigma percentile corrections
    tpcorrps = [np.percentile(tpcorrs, p, axis=0) for p in percentiles]
    opt_tpcorrps = [np.percentile(opt_tpcorrs, p, axis=0) for p in percentiles]
    
    purps = mpl.cm.ScalarMappable(mpl.colors.Normalize(vmin=0,vmax=1),'Purples')
    purps.to_rgba(1)

    # plot correction
    mpl.rcParams['font.size'] = 12
    plt.figure(figsize=(16,4))
    print wave.shape, tpcorrps[0].shape
    plt.fill_between(wave, tpcorrps[0], tpcorrps[4], color=purps.to_rgba(.25))
    plt.fill_between(wave, tpcorrps[1], tpcorrps[3], color=purps.to_rgba(.75))

    plt.plot(wave, opt_tpcorrps[0], color='black', lw=2, ls='--')
    plt.plot(wave, opt_tpcorrps[1], color='black', lw=2, ls='-')
    plt.plot(wave, opt_tpcorrps[3], color='black', lw=2, ls='-')
    plt.plot(wave, opt_tpcorrps[4], color='black', lw=2, ls='--')

    plt.axvline(4000, ls=':', c='b')
    plt.axvline(5400, ls=':', c='r')
    plt.axvline(4539, ls=':', c='k')

    plt.xlim(3500, 10500)
    # plt.grid(True)
    plt.xlabel('Observed Wavelength $\lambda$ $(\AA)$')
    plt.ylabel(r'Calibration Corrections $R_i^{\ast}(\lambda)$')
    plt.ylim(.5, 2.2)
    plt.savefig('correction_summary.pdf', bbox_inches='tight', dpi=200)
    # save figure

    design_alt_ylim = (45,90)
    delta_ha_ylim = (-35,35)
    delta_alt_lim = (-12,12)
    vmin, vmax = delta_alt_lim

    cbar_alt_ticks = np.linspace(-12,12,9)
    cbar_alt_labels = ['< -12', '-9', '-6', '-3', '0', '3', '6', '9', '> 12']

    delta_alt = mean_alt - design_alt

    mpl.rcParams['font.size'] = 18
    fig = plt.figure(figsize=(8,6))
    ax4 = plt.subplot(111)
    plt.scatter(delta_alt, delta_ha, marker='.', edgecolor='none', c=design_ha, vmin=-60, vmax=60)
    plt.ylabel('$h_\mathrm{obs} - h_\mathrm{0}$ $(\mathrm{degrees})$')
    plt.xlabel('$a_\mathrm{obs} - a_\mathrm{0}$ $(\mathrm{degrees})$')
    plt.ylim(delta_ha_ylim)
    plt.xlim(-25,25)
    plt.grid(True)
    plt.axhline(+dha_threshold, ls='--', color='k')
    plt.axhline(-dha_threshold, ls='--', color='k')

    # ax4.xaxis.set_major_locator(MultipleLocator(.25))
    plt.tick_params(axis='both', which='major', labelsize=16)

    cbar = plt.colorbar(label='$h_\mathrm{0}$ $(\mathrm{degrees})$')
    # cbar.set_ticks(cbar_alt_ticks)
    # cbar.set_ticklabels(cbar_alt_labels)
    cbar.solids.set_rasterized(True)
    cbar.solids.set_edgecolor("face")

    plt.savefig('dalt-dha-scatter.pdf')

    plt.figure(figsize=(8,6))
    plt.hist(delta_ha, bins=np.linspace(-45,45,46), alpha=0.5, histtype='stepfilled')
    plt.xlim(-45,45)
    plt.ylim(0,250)
    plt.xlabel(r'$h_{obs} - h_{design}$ (degrees)')
    plt.ylabel('Observations')
    plt.grid(True)
    plt.savefig('dha-dist.pdf')

if __name__ == '__main__':
    main()
