#!/usr/bin/python
import os
import json
import sqlite3

# THIS SCRIPT generate from 
# USAGE : python get_gene_transcript_list.py > $NGS_PIPELINE_BX_DIR/reference_files/annotate_bed/favorite_NM_list.tsv
# to be used with annotate_bed.py

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

db_cur.execute("SELECT * FROM Gene")
db_genes = db_cur.fetchall()

for db_gene in db_genes:
	print '%s\t%s' % (db_gene['geneID'],db_gene['transcript'])
