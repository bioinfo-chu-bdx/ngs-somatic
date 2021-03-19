#!/usr/bin/python
import os
import sys
import json
import sqlite3

# THIS SCRIPT IS DESIGNED TO DELETE TRANSCRIPT AND VARIANTANNOTATION ENTRIES WHEN TRANSCRIPT IS NO LONGER PRESENT IN DATABASE

# USAGE : python clean_variantbase.py

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

db_cur.execute("SELECT * FROM Transcript")
db_transcripts = db_cur.fetchall()
for db_transcript in db_transcripts:
	transcript = db_transcript['transcriptID']
	db_cur.execute("SELECT * FROM TargetedRegion WHERE transcript='%s'" % transcript)
	if db_cur.fetchone() is None:
		print "- transcript %s (%s) not longer used in DB, removing corresponding Transcript and VariantAnnotations entries" % (transcript,db_transcript['gene'])
		db_cur.execute("DELETE FROM VariantAnnotation WHERE transcript='%s'" % transcript)
		db_cur.execute("DELETE FROM Transcript WHERE transcriptID='%s'" % transcript)

# db_con.commit()
# db_con.close()
