#!/usr/bin/env python
"""
"""
import argparse
import os

import h5py
import numpy as np

import astropy.table
from astropy.io import fits

import scipy.interpolate

import bossdata

import specsim

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tpcorr', type=str, default=None,
        help='throughput correction filename, required')
    parser.add_argument('--output', type=str, default=None,
        help='output filename')
    args = parser.parse_args()

    # Open connection spAll db
    meta_db = bossdata.meta.Database(lite=False, verbose=True)

    # Find the plates used for "validation", the plates with at least 10 offset standards
    # in each spectro graph.

    # Select offset standards via ancillary target flag
    where = 'ANCILLARY_TARGET2&(1<<20)!=0'
    what = 'PLATE,MJD,FIBER,PSFMAG_1,PSFMAG_2,PSFMAG_3'
    sql = 'SELECT {} FROM meta WHERE {}'.format(what, where)
    meta_db.cursor.execute(sql)
    rows = meta_db.cursor.fetchall()
    table = astropy.table.Table(rows=rows, names=what.split(','))

    # Group by observations
    table_by_obs = table.group_by(['PLATE','MJD'])

    # Count number of targets in each spectrograph
    counts_per_spec = [(grp['PLATE'][0], np.sum(grp['FIBER'] <= 500), np.sum(grp['FIBER'] > 500)) for grp in table_by_obs.groups]
    at_least_10 = [obs[0] for obs in counts_per_spec if obs[1]+obs[2] >= 10]
    validiation_plates = [obs[0] for obs in counts_per_spec if obs[1] >= 10 and obs[2] >= 10]

    # Summarize various target sample selections used in paper
    sample_names = ['Failed quasars', 'Spec. standards', 'Offset standards']
    sample_selections = [
        'LAMBDA_EFF=4000 and OBJTYPE="QSO" and CLASS="STAR"',
        'LAMBDA_EFF=5400 and OBJTYPE="SPECTROPHOTO_STD" and CLASS="STAR"',
        'LAMBDA_EFF=4000 and CLASS="STAR" and ANCILLARY_TARGET2=(1<<20)',
    ]

    # Require targets were on "validation" plates
    validiation_plates_str = ','.join(['{}'.format(plate) for plate in validiation_plates])
    valid_selection = 'ZWARNING=0 and PLATE in ({})'.format(validiation_plates_str)
    
    # Loop over target samples
    sql_prefix = 'SELECT {} FROM meta'.format(what)
    target_lists = []
    for sample_name, sample_selection in zip(sample_names, sample_selections):
        sample_nums = {}
        
        # Count the number of targets in this sample+category
        sql = sql_prefix + ' WHERE {} and {}'.format(sample_selection, valid_selection)
        meta_db.cursor.execute(sql)
        rows = meta_db.cursor.fetchall()
        print len(rows)
        target_lists.append(rows)

    # Open the throughput correction file
    tpcorr = h5py.File(args.tpcorr, 'r')
    tpcorr_wave = tpcorr['wave'].value

    dr12_finder = bossdata.path.Finder(sas_path='/sas/dr12/boss', redux_version='v5_7_0')
    blue_finder = bossdata.path.Finder(sas_path='/sas/dr12/boss', redux_version='test')
    mirror = bossdata.remote.Manager()

    corrections = [0.036, 0.015, 0.013]

    for name, target_list in zip(sample_names, target_lists):
        for plate,mjd,fiber,psf1,psf2,psf3 in target_list:

            dr12_spec_name = dr12_finder.get_spec_path(plate, mjd, fiber=fiber, lite=True)
            dr12_spec = bossdata.spec.SpecFile(mirror.get(dr12_spec_name, timeout=None, progress_min_size=0))

            data = dr12_spec.get_valid_data(fiducial_grid=True, use_ivar=True)
            wlen,flux,ivar = data['wavelength'][:],data['flux'][:],data['ivar'][:]

            spectrum = specsim.spectrum.SpectralFluxDensity(wlen, flux)
            print psf1,psf2,psf3, spectrum.getABMagnitudes()


if __name__ == '__main__':
    main()
