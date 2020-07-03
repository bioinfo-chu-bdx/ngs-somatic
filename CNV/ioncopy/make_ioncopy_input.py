#!/usr/bin/env python
import subprocess
import sys
import os
import json
import csv
import glob
from optparse import OptionParser
import zipfile
import shutil
import sqlite3

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

############################################################################################
parser = OptionParser()
parser.add_option('-i', '--run-folder', help='Run folder ', dest='run_folder')
(options, args) = parser.parse_args()

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()
if os.path.isfile(options.run_folder+'/barcodes.json'):
	with open(options.run_folder+'/barcodes.json', 'r') as g:
		barcodes_json = json.load(g)
else:
	print "error : barcodes.json not found in run folder"

control_names = ['H2O','NTC','TEMOINCNV'] # liste des noms possibles pour les temoins negatifs
cna_dir = options.run_folder + '/_CNA'
############################################################################################

bamlist = glob.glob(options.run_folder+'/*/*.bam')
bamlist = [item for item in bamlist if not 'processed' in item]

barcode_data = {}
for bamfile in bamlist:
	if '_IonXpress' in bamfile:
		sample = bamfile.split('/')[-1].split('_IonXpress')[0]
		barcode = 'IonXpress_' + bamfile.split('IonXpress_')[-1].split('.bam')[0]
	else:
		sample = bamfile.split('/')[-1].split('_S')[0]
		barcode = 'S%s' % bamfile.split('_S')[-1].split('.bam')[0]
	barcode_data[barcode] = {}
	barcode_data[barcode]['sample'] = sample
	target = barcodes_json[barcode]['target_region_filepath'].split('/')[-1]
	for rt in global_param['run_type']:
		if global_param['run_type'][rt]['target_bed'].split('/')[-1] == target:
			barcode_data[barcode]['runtype'] = rt
			barcode_data[barcode]['targetbed'] = global_param['run_type'][rt]['target_bed']
			break
	intermediate_folder = '%s/%s/intermediate_files' % (options.run_folder, sample)
	
	if os.path.isfile('%s/%s/intermediate_files/coverage/%s_%s.amplicon.cov.xls' % (options.run_folder,sample,sample,barcode)):
		cov_file = open('%s/%s/intermediate_files/coverage/%s_%s.amplicon.cov.xls' % (options.run_folder,sample,sample,barcode),'r')
	elif os.path.isfile('%s/%s/intermediate_files/coverage/%s_%s.target.cov.xls' % (options.run_folder,sample,sample,barcode)):
		cov_file = open('%s/%s/intermediate_files/coverage/%s_%s.target.cov.xls' % (options.run_folder,sample,sample,barcode),'r')
	elif os.path.isfile('%s/%s/intermediate_files.zip' % (options.run_folder,sample)):
		archive = zipfile.ZipFile('%s/%s/intermediate_files.zip' % (options.run_folder,sample), 'r')
		if 'coverage/%s_%s.amplicon.cov.xls' % (sample,barcode) in archive.namelist():
			cov_file = archive.open('coverage/%s_%s.amplicon.cov.xls' % (sample,barcode))
		elif 'coverage/%s_%s.target.cov.xls' % (sample,barcode) in archive.namelist():
			cov_file = archive.open('coverage/%s_%s.target.cov.xls' % (sample,barcode))
	barcode_data[barcode]['covfile'] = cov_file

################################
# PROCESS POUR CHAQUE RUN TYPE #
################################

if not os.path.isdir(cna_dir):
	subprocess.call(['mkdir',cna_dir])

runtypes = list(set([barcode_data[b]['runtype'] for b in barcode_data.keys()]))
for runtype in runtypes:
	cna_runtype_dir = '%s/%s' % (cna_dir,runtype)
	if not os.path.isdir(cna_runtype_dir):
		subprocess.call(['mkdir',cna_runtype_dir])

	##############################
	# GET COVERAGE ANALYSIS DATA #
	##############################

	# TODO : remplacer "panel" par le target dans la requete ci dessous
	panel = global_param['run_type'][rt]['target_bed'].split('/')[-1]
	db_cur.execute("SELECT chromosome,start,stop,targetedRegionName,gene,details FROM TargetedRegion INNER JOIN Panel ON TargetedRegion.panel = Panel.panelID WHERE panel='%s' ORDER BY start" % panel)
	db_target_regions = db_cur.fetchall()
	region2gene = {}
	for db_target_region in db_target_regions:
		region2gene[db_target_region['targetedRegionName']] = '%s' % db_target_region['gene']

	amplicons = {}
	sample_list = []
	barcode_list = []
	for barcode in barcode_data.keys():
		sample = barcode_data[barcode]['sample']
		if barcode_data[barcode]['runtype'] != runtype: # check if barcode is in runtype
			continue
		iscontrol = False # check if barcode is not a control
		for control_name in control_names:
			if control_name in sample.upper():
				iscontrol = True
		if iscontrol:
			continue
		sample_list.append(sample)
		barcode_list.append(barcode)
		# ouverture du fichier de couverture par amplicon de chaque patient, recuperation des infos
		reader = csv.reader(barcode_data[barcode]['covfile'], delimiter = '\t')
		reader.next()
		for row in reader:
			amplicon_id = row[3]
			gene = region2gene[amplicon_id]
			amplgene = gene + '_' + amplicon_id
			if not amplgene in amplicons:
				amplicons[amplgene] = {}
			total_reads = row[9]
			amplicons[amplgene][barcode] = total_reads
		intermediate_folder = '%s/%s/intermediate_files' % (options.run_folder, sample)
		# re-zipper uniquement si il a fallu dezippe
		#shutil.make_archive(intermediate_folder,'zip',intermediate_folder)
		#shutil.rmtree(intermediate_folder)

	#########################
	## WRITE IONCOPY INPUT ##
	#########################

	incopy_input_file = open('%s/ioncopy_input.tsv' % cna_runtype_dir,'w')
	incopy_input_writer = csv.writer(incopy_input_file,delimiter='\t')
	
	incopy_input_writer.writerow([''] + sample_list)
	
	for amplgene in amplicons.keys():
		line = [amplgene]
		for barcode in barcode_list:
			# IONCOPY bug sur valeurs nombre de reads = 0. Remplacer les 0 par des 1 ne change pas le resultat et permet d'eviter le bug
			if str(amplicons[amplgene][barcode]) == '0':
				line.append('1')
			else:
				line.append(amplicons[amplgene][barcode])
		incopy_input_writer.writerow(line)
	incopy_input_file.close()
