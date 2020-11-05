#!/usr/bin/python
import os
import sys
import csv
import json
import sqlite3
import subprocess
import uuid
import hgvs.dataproviders.uta

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

with open(global_param['NC_data'], 'r') as ncdata:
	NC_data = json.load(ncdata)

refGene_file = open(global_param['RefSeq'],'r')
refGene_reader = csv.reader(refGene_file,delimiter='\t')
hdp = hgvs.dataproviders.uta.connect()
db_path = global_param['VariantBase']

db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

refGene_file.seek(0)

for transcript in ['NM_000516.5','NM_001354870.1','NM_001025203.1','NM_001198551.1','NM_022552.4']:
for rgline in refGene_reader:
	if rgline[1].split('.')[0] == transcript and ('_' not in rgline[2]) :
		transcriptionStart = rgline[4]
		transcriptionStop = rgline[5]	# rgline[6] et rgline[7] sont CodingRegionStart et codingRegionStop
		exons = int(rgline[8])
		exonsStart = rgline[9]
		exonsStop =  rgline[10]
		break

gene_tx = hdp.get_tx_for_gene(gene)
version = 0
for item in gene_tx:
	if (item[3].split('.')[0] == transcript) and (item[4] == nc):
		version = max(version,int(item[3].split('.')[-1]))	
db_cur.execute("INSERT INTO Gene (geneID, chromosome, transcriptionStart, transcriptionStop, strand, exons, exonsStart, exonsStop, transcript, transcriptVersion, NC) VALUES ('%s', '%s', %s, %s, '%s', %s, '%s', '%s', '%s', %s, '%s')" % (gene,chromosome,transcriptionStart,transcriptionStop,strand,exons,exonsStart,exonsStop,transcript,version,nc))

