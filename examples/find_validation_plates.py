#!/usr/bin/env python

import numpy as np

import astropy.table

import bossdata.meta

def main():

	meta_db = bossdata.meta.Database(lite=False, verbose=True)
	meta_db.cursor.execute('SELECT PLATE,MJD,FIBER FROM meta WHERE ((ANCILLARY_TARGET2&(1<<20))>0)')
	rows = meta_db.cursor.fetchall()
	table = astropy.table.Table(rows=rows, names=['PLATE','MJD','FIBER'])
	table_by_obs = table.group_by(['PLATE','MJD'])

	counts_per_spec = [(grp['PLATE'][0], np.sum(grp['FIBER'] <= 500),np.sum(grp['FIBER'] > 500)) for grp in table_by_obs.groups]
	at_least_10 = [obs[0] for obs in counts_per_spec if obs[1]+obs[2] >= 10]
	validiation_plates = [obs[0] for obs in counts_per_spec if obs[1] >= 10 and obs[2] >= 10]

	print 'Number of observations with ancillary targets:', len(counts_per_spec)
	print 'Number of observations with at least 10 ancillary targets: ', len(at_least_10)
	print 'Number of observations with at least 10 ancillary targets per spectrograph:', len(validiation_plates)
	print 'Validation plates:', validiation_plates

	bad_chunks = ('boss35','boss36','boss37','boss38')

	# LAMBDA_EFF=4000 and ZWARNING&1<<7=0 and CHUNK not in ('boss35','boss36','boss37','boss38')

	validiation_plates_str = ','.join(['{}'.format(plate) for palte in validiation_plates])

	meta_db.cursor.execute('SELECT PLATE,MJD,FIBER FROM meta WHERE LAMBDA_EFF=4000 and ZWARNING=0 and PLATE in ({})'.format(
		validiation_plates_str))
	rows = meta_db.cursor.fetchall()
	print 'Number of observations with offset target on validation plates:', len(rows)

if __name__ == '__main__':
	main()