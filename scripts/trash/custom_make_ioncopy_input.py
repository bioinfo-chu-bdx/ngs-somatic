#!/usr/bin/env python
import subprocess
import sys
import os
import json
import csv
import glob
from optparse import OptionParser

############################################################################################
parser = OptionParser()
parser.add_option('-i', '--run-folder', help='Run folder ', dest='run_folder')
(options, args) = parser.parse_args()

with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)

if os.path.isfile(options.run_folder+'/barcodes.json'):
	with open(options.run_folder+'/barcodes.json', 'r') as g:
		barcodes_json = json.load(g)
else:
	print "error : barcodes.json not found in run folder"

control_names = ['H2O','H20','NTC','TEMOIN NEG','EAU','TEMOINCNV'] # liste des noms possibles pour les temoins negatifs
cna_dir = options.run_folder + '/_2018_reanalysis/_CNA'
############################################################################################

bamlist = glob.glob(options.run_folder+'/*/*.bam')
bamlist = [item for item in bamlist if not 'processed' in item]

barcode2runtype = {}
barcode2sample = {}
barcode2covfile = {}
for bamfile in bamlist:
	sample = bamfile.split('/')[-1].split('_IonXpress')[0]
	barcode = 'IonXpress_' + bamfile.split('IonXpress_')[-1].split('.bam')[0]
	barcode2sample[barcode] = sample
	target = barcodes_json[barcode]['target_region_filepath'].split('/')[-1]
	for _run_type in global_param['run_type']:
		if global_param['run_type'][_run_type]['target_bed'].split('/')[-1] == target:
			barcode2runtype[barcode] = _run_type
			break
	covfile = '/DATA/work/results/%s_%s.bam/coverage/%s_%s.amplicon.cov.xls' % (sample.replace(' ','@'),barcode,sample,barcode)
	if os.path.isfile(covfile):
		barcode2covfile[barcode] = covfile

################################
# PROCESS POUR CHAQUE RUN TYPE #
################################

if not os.path.isdir(cna_dir):
	subprocess.call(['mkdir',cna_dir])

for runtype in list(set(barcode2runtype.values())):
	cna_runtype_dir = '%s/%s' % (cna_dir,runtype)
	if not os.path.isdir(cna_runtype_dir):
		subprocess.call(['mkdir',cna_runtype_dir])
		
	##############################
	# GET COVERAGE ANALYSIS DATA #
	##############################

	amplicons = {}
	sample_list = []
	barcode_list = []
	for barcode in barcode2covfile.keys():
		sample = barcode2sample[barcode]
		if barcode2runtype[barcode] != runtype: # check if barcode is in runtype
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
		with open(barcode2covfile[barcode],'r') as cov_file:
			reader = csv.reader(cov_file, delimiter = '\t')
			reader.next()
			for row in reader:
				amplicon_id = row[3]
				gene_id = row[4].split('GENE_ID=')[-1].split(';')[0]
				amplgene = gene_id + '_' + amplicon_id
				if not amplgene in amplicons:
					amplicons[amplgene] = {}
				total_reads = row[9]
				amplicons[amplgene][barcode] = total_reads

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
