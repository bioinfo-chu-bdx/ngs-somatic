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

# DELETE EVERYTHING RELATED TO A RUN : ANALYSIS, VARIANTMETRICS AND RUN.
run = 'Auto_user_S5-0198-470-NGS518_519_35pM_Chef_SBT-colon-lung_v10_530_700'
db_cur.execute("DELETE FROM VariantMetrics WHERE variantMetricsID IN (SELECT variantMetricsID FROM VariantMetrics INNER JOIN Analysis ON Analysis.analysisID = VariantMetrics.analysis INNER JOIN Run ON Run.runID = Analysis.run WHERE Run.runID='%s')" % run)
db_cur.execute("DELETE FROM Analysis WHERE analysisID IN (SELECT analysisID FROM Analysis INNER JOIN Run ON Run.runID = Analysis.run WHERE Run.runID='%s')" % run)
db_cur.execute("DELETE FROM Run WHERE Run.runID='%s'" % run)

# db_cur.execute("SELECT * FROM VariantMetrics")
# db_trucs = db_cur.fetchall()
# for db_truc in db_trucs:
	# corrected_truc = False
	# truc_id = db_truc['variantMetricsID']
	# call = db_truc['call']
	# if db_truc['variantCallingTool'] == 'TVC':
		# corrected_truc = call.replace(' / ','/')
		# corrected_truc = corrected_truc.replace('de novo','tvc_de_novo')
		# corrected_truc = corrected_truc.replace('hotspot','tvc_hotspot')
	# if corrected_truc:
		# db_cur.execute("UPDATE VariantMetrics SET call='%s' WHERE variantMetricsID='%s'" % (corrected_truc,truc_id))

# db_cur.execute("SELECT * FROM VariantAnnotation")
# db_trucs = db_cur.fetchall()
# for db_truc in db_trucs:
	# truc_id = db_truc['variantAnnotationID']
	# truc_date = db_truc['lastUpdate']
	# if '/' in truc_date:
		# truc_date = truc_date.split('/')
		# j = truc_date[0]
		# m = truc_date[1]
		# a = truc_date[2]
	# elif '-' in truc_date:
		# truc_date = truc_date.split('-')
		# j = truc_date[2]
		# m = truc_date[1]
		# a = truc_date[0]
	# corrected_truc_date = '%s%s%s' % (a,m,j)
	# db_cur.execute("UPDATE VariantAnnotation SET lastUpdate='%s' WHERE variantAnnotationID='%s'" % (corrected_truc_date,truc_id))

db_con.commit()
db_con.close()
