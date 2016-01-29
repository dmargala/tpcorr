#!/usr/bin/env python

import numpy as np

import astropy.table

import bossdata.meta

def main():

    # Open connection spAll db
    meta_db = bossdata.meta.Database(lite=False, verbose=True)

    # Find the plates used for "validation", the plates with at least 10 offset standards
    # in each spectro graph.

    # Select offset standards via ancillary target flag
    where = 'ANCILLARY_TARGET2&(1<<20)!=0'
    what = 'PLATE,MJD,FIBER'
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

    # Print summary
    print 'Number of observations with ancillary targets: {}'.format(len(counts_per_spec))
    print 'Number of observations with at least 10 ancillary targets: {}'.format(len(at_least_10))
    print 'Number of observations with at least 10 ancillary targets per spectrograph: {}'.format(len(validiation_plates))
    print 'Validation plates: {}'.format(validiation_plates)

    # Summarize various target sample selections used in paper
    sample_names = ['Offset targets', 'Quasars', 'Failed quasars', 'Spec. standards', 'Offset standards']
    sample_selections = [
        'LAMBDA_EFF=4000',
        'LAMBDA_EFF=4000 and OBJTYPE="QSO" and CLASS="QSO"',
        'LAMBDA_EFF=4000 and OBJTYPE="QSO" and CLASS="STAR"',
        'LAMBDA_EFF=5400 and OBJTYPE="SPECTROPHOTO_STD" and CLASS="STAR"',
        'LAMBDA_EFF=4000 and CLASS="STAR" and ANCILLARY_TARGET2=(1<<20)',
    ]

    print '\nTarget sample definitions: '
    for sample_name, sample_selection in zip(sample_names, sample_selections):
        print '{}: {}'.format(sample_name, sample_selection)

    # Construct category selections
    category_names = ('DR12', 'DR12b', 'Validation', 'DR12c')

    # Require that the fiber is plugged
    dr12_selection = 'ZWARNING&(1<<7)=0'

    # "Bad chunk" observations have spectrophoto standards with focal plane offsets. Offset targets
    # on these plates therefore do not require a correction
    bad_chunks = ('boss35','boss36','boss37','boss38')
    bad_chunks_str = ','.join(['"{}"'.format(chunk) for chunk in bad_chunks])
    dr12b_selection = 'ZWARNING&(1<<7)=0 and CHUNK not in ({})'.format(bad_chunks_str)
    
    dr12c_selection = 'ZWARNING=0 and CHUNK in ({})'.format(bad_chunks_str)

    # Require targets were on "validation" plates
    validiation_plates_str = ','.join(['{}'.format(plate) for plate in validiation_plates])
    valid_selection = 'ZWARNING=0 and PLATE in ({})'.format(validiation_plates_str)

    # Print category definitions
    category_selections = (dr12_selection, dr12b_selection, valid_selection, dr12c_selection)
    print '\nCategory definitions: '
    for category_name, category_selection in zip(category_names, category_selections):
        print '{}: {}'.format(category_name, category_selection)

    # Loop over target samples
    print '\nCounts: '
    sql_prefix = 'SELECT {} FROM meta'.format(what)
    category_nums = {}
    for sample_name, sample_selection in zip(sample_names, sample_selections):
        sample_nums = {}
        # Loop over sample categories
        for category_name, category_selection in zip(category_names, category_selections):
            # Count the number of targets in this sample+category
            sql = sql_prefix + ' WHERE {} and {}'.format(sample_selection, category_selection)
            meta_db.cursor.execute(sql)
            rows = meta_db.cursor.fetchall()
            sample_nums[category_name] = len(rows)
            # For the DR12 offset target sample, save number of obs and plates
            if sample_name == sample_names[0]:
                table = astropy.table.Table(rows=rows, names=what.split(','))
                table_by_obs = table.group_by(['PLATE','MJD'])
                table_by_plate = table.group_by(['PLATE'])
                category_nums[category_name] = (len(table_by_obs.groups), len(table_by_plate.groups))
        # Print sample summary
        print '{}: {}'.format(sample_name, ', '.join(
            [str(sample_nums[category_name]) for category_name in category_names]))

    # Print number of obs and plates with offset targets
    print '\nNumber observation (unique plate) with offset targets summary:'
    for category_name in category_names:
        print '{}: {} ({})'.format(category_name, category_nums[category_name][0], category_nums[category_name][1])

if __name__ == '__main__':
    main()
