#!/usr/bin/python
import subprocess
import sqlite3
import pysam
import glob
import json
import csv
import re
import os

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

# liste des noms possibles pour les temoins negatifs
control_names = ['H2O','H20','NTC','ACROMETRIX','BAF5','BAF-5','HD300','HD301','HD802','HD901','HD748','HORIZON','TEMOIN']
sbt_pattern = r"^.*[A-Z]{2}[0-9]{3}[A-Z]{1}$" ## ex BESCOND-AW210F (2 lettres 3 chiffres 1 lettre)
sbt_old_pattern = r"^.*[A-Z]{1}[0-9]{3}[A-Z]{1}$" ## ex GUSMINI-Z492F (1 lettre 3 chiffres 1 lettre)

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
db_con = sqlite3.connect('%s/variantBase/VariantBase.db' % pipeline_folder)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()


# db_cur.execute("SELECT * FROM Run WHERE platform='ion torrent'")
# db_cur.execute("SELECT * FROM Run WHERE platform='ion torrent' AND runID='Auto_user_S5-0198-470-NGS518_519_35pM_Chef_SBT-colon-lung_v10_530_700'")
# db_trucs = db_cur.fetchall()

db_trucs = [{'runID':'Auto_user_S5-0198-497-NGS534_535_35pM_Chef_SBT-colon-lung_v10_530_727','runPath':'/media/n06lbth/sauvegardes_pgm/SBT/Run_500-599/Auto_user_S5-0198-497-NGS534_535_35pM_Chef_SBT-colon-lung_v10_530_727','system':'S5'},{'runID':'Auto_user_S5-0198-498-Run223-TP53-ABL-FLT3-MM_729','runPath':'/media/n06lbth/sauvegardes_pgm/LAM/panel-myeloid-v8/Auto_user_S5-0198-498-Run223-TP53-ABL-FLT3-MM_729','system':'S5'}]

for db_truc in db_trucs:
	barcodes_json_path = '%s/barcodes.json' % db_truc['runPath']
	if not os.path.exists(barcodes_json_path):
		print "barcodes.json not found for run %s" % db_truc['runID']
		continue
	print "- RUN : %s" % db_truc['runID']
	with open(barcodes_json_path, 'r') as g:
		barcodes_json = json.load(g)

	new_barcodes_json = {}

	for barcode in barcodes_json:
		# KEEP
		new_barcodes_json[barcode] = {}
		new_barcodes_json[barcode]['sample'] = barcodes_json[barcode]['sample']
		new_barcodes_json[barcode]['sample_id'] = barcodes_json[barcode]['sample_id']
		new_barcodes_json[barcode]['reference'] = 'hg19'

		# TRANSFORM
		new_barcodes_json[barcode]['library'] = barcodes_json[barcode]['barcode_description']
		new_barcodes_json[barcode]['target_bed'] = barcodes_json[barcode]['target_region_filepath'].split('/')[-1]
		
		# ADD
		new_barcodes_json[barcode]['analysis_id'] = ''
		new_barcodes_json[barcode]['description'] = ''
		new_barcodes_json[barcode]['platform'] = 'ion torrent'
		new_barcodes_json[barcode]['system'] = db_truc['system']
		new_barcodes_json[barcode]['target_technique'] = 'amplicon'

		is_control = 0
		for control_name in control_names:
			if control_name in new_barcodes_json[barcode]['sample'].upper():
				# control should not have dna number... this test is to avoid case like : REVLEG-AH200F
				if re.match(sbt_pattern,new_barcodes_json[barcode]['sample']) or re.match(sbt_old_pattern,new_barcodes_json[barcode]['sample']):
					continue
				else:
					is_control = 1
		new_barcodes_json[barcode]['is_control'] = is_control

		new_barcodes_json[barcode]['panel'] = ''
		if 'IAD172906' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'ColonLung_v10'
		elif 'IAD165023' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'ColonLung_v9'
		elif 'IAD154118' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'ColonLung_v8'
		elif 'IAD119108' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'ColonLung_v7'
		elif 'IAD108862' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'ColonLung_v5'
		elif 'IAD94971' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'ColonLung_v4'
		elif 'IAD72953' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'ColonLung_v3'
		elif 'TP53' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'TP53'
		elif 'IAD119887' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'Lymphome_B'
		elif 'IAD120574' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'Lymphome_T'
		elif 'IAD83112' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'Leuc'
		elif 'ABL1' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'ABL1'
		elif 'FLT3' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'FLT3'
		elif 'IAD78219' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'LAM-iontorrent-v1'
		elif 'IAD37093' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'LAM-iontorrent-v2'
		elif 'IAD62716' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'LAM-iontorrent-v3'
		elif 'IAD143291' in barcodes_json[barcode]['target_region_filepath']:
			new_barcodes_json[barcode]['panel'] = 'LAM-iontorrent-v8'

	# BACKUP OLD BARCODES JSON
	subprocess.call(['mv',barcodes_json_path,'%s/barcodes.old.json' % db_truc['runPath']])

	print "- writing new barcodes JSON..."
	json_text = json.dumps(new_barcodes_json, indent=4, sort_keys=True)
	bc_json = open(barcodes_json_path,'w')
	bc_json.write(json_text)
	bc_json.close()
