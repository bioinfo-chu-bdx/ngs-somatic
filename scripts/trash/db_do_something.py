#!/usr/bin/python
import sqlite3
import pysam
import json
import csv
import os

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
db_con = sqlite3.connect('%s/variantBase/VariantBase.db' % pipeline_folder)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()


db_cur.execute("SELECT * FROM Analysis WHERE panel in ('LAM-illumina-v1','LAM-illumina-v2')")
db_trucs = db_cur.fetchall()
for db_truc in db_trucs:
	bamPath = db_truc['bamPath']
	if 'LAM/panel-capture-v1' in bamPath :
		bamPath = bamPath.replace('LAM/panel-capture-v1','LAM/panel-capture')
	if 'Run1-Test-SureSelect-Myeloid' in bamPath :
		bamPath = bamPath.replace('Run1-Test-SureSelect-Myeloid','Run-01-Test-SureSelect-Myeloid')
	if 'Run2-Test-SureSelect-Myeloid' in bamPath :
		bamPath = bamPath.replace('Run2-Test-SureSelect-Myeloid','Run-02-Test-SureSelect-Myeloid')
	if 'Run-3-CAP-MYELOID' in bamPath :
		bamPath = bamPath.replace('Run-3-CAP-MYELOID','Run-03-CAP-MYELOID')
	if 'Run-5-CAP-MYELOID' in bamPath :
		bamPath = bamPath.replace('Run-5-CAP-MYELOID','Run-05-CAP-MYELOID')
	db_cur.execute("UPDATE Analysis SET bamPath='%s' WHERE analysisID='%s'" % (bamPath,db_truc['analysisID']))

db_con.commit()
db_con.close()
