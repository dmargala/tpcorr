#!/usr/bin/env python

import numpy as np

import astropy.table

import bossdata.meta

def main():

    meta_db = bossdata.meta.Database(lite=False, verbose=True)

    what = 'PLATE,MJD,FIBER'
    where = 'ANCILLARY_TARGET2&(1<<20)!=0'
    sql = 'SELECT {} FROM meta WHERE {}'.format(what, where)
    meta_db.cursor.execute(sql)
    rows = meta_db.cursor.fetchall()
    table = astropy.table.Table(rows=rows, names=what.split(','))
    table_by_obs = table.group_by(['PLATE','MJD'])

    counts_per_spec = [(grp['PLATE'][0], np.sum(grp['FIBER'] <= 500),np.sum(grp['FIBER'] > 500)) for grp in table_by_obs.groups]
    at_least_10 = [obs[0] for obs in counts_per_spec if obs[1]+obs[2] >= 10]
    validiation_plates = [obs[0] for obs in counts_per_spec if obs[1] >= 10 and obs[2] >= 10]

    print 'Number of observations with ancillary targets:', len(counts_per_spec)
    print 'Number of observations with at least 10 ancillary targets: ', len(at_least_10)
    print 'Number of observations with at least 10 ancillary targets per spectrograph:', len(validiation_plates)
    print 'Validation plates:', validiation_plates

    bad_chunks = ('boss35','boss36','boss37','boss38')
    bad_chunks_str = ','.join(['"{}"'.format(chunk) for chunk in bad_chunks])

    # LAMBDA_EFF=4000 and ZWARNING&1<<7=0 and CHUNK not in ('boss35','boss36','boss37','boss38')

    validiation_plates_str = ','.join(['{}'.format(plate) for plate in validiation_plates])

    sample_names = ['Offset targets', 'Quasars', 'Failed quasars', 'Spec. standards', 'Offset standards']
    sample_selections = [
        'LAMBDA_EFF=4000 and ZWARNING&(1<<7)=0',
        'LAMBDA_EFF=4000 and ZWARNING&(1<<7)=0 and OBJTYPE="QSO" and CLASS="QSO"',
        'LAMBDA_EFF=4000 and ZWARNING&(1<<7)=0 and OBJTYPE="QSO" and CLASS="STAR"',
        'LAMBDA_EFF=5400 and ZWARNING&(1<<7)=0 and OBJTYPE="SPECTROPHOTO_STD" and CLASS="STAR"',
        'LAMBDA_EFF=4000 and ZWARNING&(1<<7)=0 and CLASS="STAR" and ANCILLARY_TARGET2=(1<<20)',
    ]

    print '\nTarget sample definitions: '

    for name, selection in zip(sample_names, sample_selections):
        print '{}: {}'.format(name, selection)

    dr12_sql_postfix = 'ZWARNING&(1<<7)=0'
    dr12b_sql_postfix = 'CHUNK not in ({})'.format(bad_chunks_str)
    valid_sql_postfix = 'ZWARNING=0 and PLATE in ({})'.format(validiation_plates_str)

    print '\nAdditional definitions: '
    print 'DR12: {}'.format(dr12_sql_postfix)
    print 'DR12b: {}'.format(dr12b_sql_postfix)
    print 'Validation: {}'.format(valid_sql_postfix)

    print '\nCounts: '
    for name, selection in zip(sample_names, sample_selections):
        sql_prefix = 'SELECT {} FROM meta'.format(what)
        dr12_sql = sql_prefix + ' WHERE {} and {}'.format(selection, dr12_sql_postfix)
        dr12b_sql = sql_prefix + ' WHERE {} and {}'.format(selection, dr12b_sql_postfix)
        valid_sql = sql_prefix + ' WHERE {} and {}'.format(selection, valid_sql_postfix)

        meta_db.cursor.execute(dr12_sql)
        rows = meta_db.cursor.fetchall()
        num_dr12 = len(rows)

        meta_db.cursor.execute(dr12b_sql)
        rows = meta_db.cursor.fetchall()
        num_dr12b = len(rows)

        meta_db.cursor.execute(valid_sql)
        rows = meta_db.cursor.fetchall()
        num_valid = len(rows)

        print '{}: {}, {}, {}'.format(name, num_dr12, num_dr12b, num_valid)

if __name__ == '__main__':
    main()
