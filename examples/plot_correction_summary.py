#!/usr/bin/env python


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
    tpcorr = h5py.File('corrections.hdf5')

    #wave = np.squeeze(tpcorr['wave'].value)
    wave = np.linspace(3500., 10500., 15)

    tpcorrs = []
    ha = []
    alt = []
    design_ha = []
    design_alt = []
    plate_mjd_list = []
    for plate in tpcorr.keys():
        if plate == 'wave':
            continue
        for mjd in tpcorr[plate].keys():
            plate_mjd_list.append((plate, mjd))
            grp = tpcorr['{}/{}'.format(plate, mjd)]
            ha.append(grp.attrs['ha'])
            design_ha.append(grp.attrs['design_ha'])
            alt.append(grp.attrs['alt'])
            design_alt.append(grp.attrs['design_alt'])
            for fiber in grp.keys():
                tpcorrs.append(tpcorr[plate][mjd][fiber].value.tolist())
    tpcorrArray = np.array(tpcorrs)

    mean_ha = np.array(ha)
    design_ha = np.array(design_ha)
    mean_alt = np.array(alt)
    design_alt = np.array(design_alt)

    dha_threshold = 1.25
    delta_ha = mean_ha - design_ha
    optimal_obs = np.where(np.abs(delta_ha) < dha_threshold)[0]
    other_obs = np.where(np.abs(delta_ha) >= dha_threshold)[0]
    # optimal_plate_mjds = zip(getField('plate')[optimal_obs].tolist(),getField('mjd')[optimal_obs])
    optimal_plate_mjd_list = [plate_mjd for plate_mjd,is_optimal in zip(plate_mjd_list, optimal_obs) if is_optimal]
    print float(len(mean_ha)-len(optimal_obs))/len(mean_ha)

    opt_tpcorrs = []
    for plate, mjd in optimal_plate_mjd_list:
        for fiber in tpcorr[plate][mjd].keys():
            opt_tpcorrs.append(tpcorr[plate][mjd][fiber].value.tolist())
    opt_tpcorrArray = np.array(opt_tpcorrs)

    # calculate cumulative dist values for a normal distribution
    percentiles = [100*scipy.stats.norm.cdf(sigma) for sigma in [-2,-1,0,1,2]]

    # Mean correction
    tpcorrmean = np.mean(tpcorrArray, axis=0)
    # N-sigma percentile corrections
    tpcorrps = [np.percentile(tpcorrArray, percentile, axis=0) for percentile in percentiles]

    # opt_tpcorrmean = np.mean(opt_tpcorrArray, axis=0)
    # opt_tpcorrps = [np.percentile(opt_tpcorrArray, percentile, axis=0) for percentile in percentiles]
    opt_tpcorrmean = np.mean(tpcorrArray, axis=0)
    opt_tpcorrps = [np.percentile(tpcorrArray, percentile, axis=0) for percentile in percentiles]
    #  5 =>  4000
    # 25 =>  6000
    # 40 =>  7500
    # -6 => 10000
    # waveIndex = 25
    waveIndex = 8
    print 'Choosing percentile examples at wavelength: ', wave[waveIndex]
    def getPercentileIndex(array, percentile):
        '''
        Returns the index of the element matching the specified percentile in array
        '''
        return np.argsort(array)[round((len(array)-1)*(percentile/100.))]

    # tpcorrmean = np.mean(tpcorrArray, axis=0)
    # tpcorrps = [np.percentile(tpcorrArray, p, axis=0) for p in sigmas]

    tpcorrIndices = [getPercentileIndex(tpcorrArray[:,waveIndex], percentile) for percentile in percentiles]
    # opt_tpcorrIndices = [getPercentileIndex(opt_tpcorrArray[:,25], percentile) for percentile in percentiles]
    purps = mpl.cm.ScalarMappable(mpl.colors.Normalize(vmin=0,vmax=1),'Purples')
    purps.to_rgba(1)

    # plot correction
    mpl.rcParams['font.size'] = 12
    plt.figure(figsize=(16,4))
    # plt.plot(wave, tpcorrmean, color='black', ls='--')
    print wave.shape, tpcorrps[0].shape
    plt.fill_between(wave, tpcorrps[0], tpcorrps[4], color=purps.to_rgba(.25))
    plt.fill_between(wave, tpcorrps[1], tpcorrps[3], color=purps.to_rgba(.75))

    plt.plot(wave, opt_tpcorrps[0], color='black', lw=2, ls='--')
    plt.plot(wave, opt_tpcorrps[1], color='black', lw=2, ls='-')
    plt.plot(wave, opt_tpcorrps[3], color='black', lw=2, ls='-')
    plt.plot(wave, opt_tpcorrps[4], color='black', lw=2, ls='--')

    # plt.plot(wave, tpcorrArray[tpcorrIndices[0]], color='blue')
    # plt.plot(wave, tpcorrArray[tpcorrIndices[2]], color='blue')
    # plt.plot(wave, tpcorrArray[tpcorrIndices[4]], color='blue')

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

    #plt.savefig('dalt-dha-scatter.pdf')

if __name__ == '__main__':
    main()
