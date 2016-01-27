#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
mpl.rcParams.update({'font.size': 18})
mpl.rcParams.update({'savefig.dpi': 200})
mpl.rcParams.update({'savefig.bbox': 'tight'})

import matplotlib.pyplot as plt

# from astropy.utils.data import download_file
# from astropy.utils import iers
# iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL, cache=True))

import astropy.units as u
import astropy.coordinates
import astropy.time

import bossdata

import os
import pydl

# helper function for shifting angles to desired range
def normalize_angle(angle):
    while angle <= -180:
        angle += 360
    while angle > 180:
        angle -= 360
    return angle

def main():

    platelist = bossdata.meta.Database(platelist=True)
    # Select columns from "good" plates
    good_plates = platelist.select_all(what='PLATE,MJD,MAPNAME,RACEN,DECCEN,SEEING50,RMSOFF50,NGUIDE,TAI_BEG,AIRTEMP',where='PLATEQUALITY="good"')

    # Some observations near the beginning of the survey do not have guide cam data
    num_missing_guide_data = np.sum(good_plates['SEEING50'] == 0)
    print 'Number of plates with SEEING50=0: %d' % num_missing_guide_data 

    # Plot seeing distribution
    plt.figure(figsize=(8,6))
    plt.hist(good_plates['SEEING50'], bins=np.linspace(0.96, 2.6, 42), alpha=.5, histtype='stepfilled')
    plt.xlim(0.96, 2.6)
    plt.ylim(0, 200)
    plt.ylabel('Observations')
    plt.xlabel('PSF FWHM (arcseconds)')
    plt.grid(True)
    plt.savefig('psf-dist.pdf')

    # Some observations near the beginning of the survey do not have guide cam data
    num_missing_rms = np.sum((good_plates['RMSOFF50'] == 0) & (good_plates['NGUIDE'] == 0))
    print 'Number of plates with RMSOFF50=0: %d' % num_missing_rms 

    # Plot seeing distribution
    plt.figure(figsize=(8,6))
    plt.hist(good_plates['RMSOFF50']*np.sqrt(good_plates['NGUIDE']), bins=np.linspace(0, 0.2, 101), alpha=.5, histtype='stepfilled')
    plt.xlim(0, 0.2)
    # plt.ylim(0, 200)
    plt.ylabel('Observations')
    plt.xlabel('RMSOFF50 (arcseconds)')
    plt.grid(True)
    plt.savefig('rms-dist.pdf')

    # It seems strange that so many plates have airtemp == 0, not sure why
    num_airtemp_zero = np.sum(good_plates['AIRTEMP'] == 0)
    print 'Number of plates with AIRTEMP=0: %d' % num_airtemp_zero 

    # Plot airtemp distribution
    plt.figure(figsize=(8,6))
    plt.hist(good_plates['AIRTEMP'], bins=np.linspace(-15, 25, 80+1), alpha=.5, histtype='stepfilled')
    plt.xlim(-15, 25)
    plt.ylabel('Observations')
    plt.xlabel('AIRTEMP (Degree Celsius)')
    plt.grid(True)
    plt.savefig('airtemp-dist.pdf')

    # observatory 
    apo = astropy.coordinates.EarthLocation.of_site('apo')
    
    tai_mid = good_plates['TAI_BEG']
    plate_centers = good_plates['RACEN']*u.deg
    when = astropy.time.Time(tai_mid/86400., format='mjd', scale='tai', location=apo)

    ha_array = np.array(map(normalize_angle, ((when.sidereal_time('apparent') - plate_centers).to(u.deg).value)))

    # finder = bossdata.path.Finder()
    # mirror = bossdata.remote.Manager()

    # speclog_path = os.getenv('BOSS_SPECLOG')

    # design_ha_array = []
    # for i,obs in enumerate(good_plates):
    #     plate,mjd = obs['PLATE'], obs['MJD']
    #     if i and (i % 25) == 0:
    #         print plate, mjd
    #     obs_mjd = '{:d}'.format(mjd)
    #     plug_map_name = 'plPlugMapM-{}.par'.format(obs['MAPNAME'])
    #     plug_map_path = os.path.join(speclog_path, obs_mjd, plug_map_name)
        
    #     try:
    #         plug_map = pydl.pydlutils.yanny.yanny(plug_map_path, np=True)
    #         # Get the list of exposures used in this observation's coadd from a spec lite file.
    #     #     spec_name = finder.get_spec_path(plate, mjd, fiber=1, lite=True)
    #     #     spec_file = bossdata.spec.SpecFile(mirror.get(spec_name))
    #         # Read the first b1 raw science exposure to find this plate's plug map.
    #     #     raw = spec_file.get_raw_image(0, 'blue', finder=finder, mirror=mirror)
    #     #     plug_map = raw.read_plug_map()

    #         # Look up the plate design pointing from the raw header and convert to
    #         # an index A,B,C,... -> 0,1,2,...
    #         # pointing_label = raw.header['POINTING'].strip()
    #         # pointing_index = ord(pointing_label) - ord('A')
    #         pointing_index = 0

    #         # Initialize a pointing object for this plate's sky location.
    #         ra_center = float(plug_map['raCen']) * u.deg
    #         dec_center = float(plug_map['decCen']) * u.deg
    #     #     print 'Plate center is RA={:.3f}, DEC={:.3f} for {}-{}'.format(ra_center, dec_center, plate, pointing_label)
    #         #pointing = tpcorr.pointing.Pointing(ra_center, dec_center)

    #         design_ha = float(plug_map['ha'].split()[pointing_index])
    #         design_ha_array.append(design_ha)
    #     except KeyError, e:
    #         print e
    #         print plate, mjd, plug_map_path

    # design_ha_array = np.array(design_ha_array)

    # plt.figure(figsize=(8,6))
    # plt.hist(ha_array-design_ha_array, bins=np.linspace(-45,45,91), alpha=0.5, histtype='stepfilled')
    # plt.xlabel(r'$h_{obs} - h_{design}$ (degrees)')
    # plt.ylabel('Observations')
    # plt.grid(True)
    # plt.savefig('dha-dist2.pdf')

if __name__ == '__main__':
    main()