#!/usr/bin/env python
import os
import sys
import json
import csv
import subprocess
import glob
import numpy
import matplotlib
matplotlib.use('Agg')
import pysam
from pylab import *
from optparse import OptionParser

"""Script searching ASXL1 dupG mutation"""

############################################################################################
parser = OptionParser()
parser.add_option('-i', '--run-folder', help='Run folder ', dest='run_folder')
(options, args) = parser.parse_args()

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

if os.path.isfile(options.run_folder+'/barcodes.json'):
	with open(options.run_folder+'/barcodes.json', 'r') as g:
		barcodes_json = json.load(g)
else:
	print "error : barcodes.json not found in run folder"

control_names = ['H2O','H20','NTC'] # liste des noms possibles pour les temoins negatifs
dupg_dir = options.run_folder + '/_dupG'

############################################################################################

bamlist = glob.glob(options.run_folder+'/*/*.bam')
bamlist = [item for item in bamlist if not 'processed' in item]

barcode2sample = {}
barcode2runtype = {}
barcode2bam = {}
for bamfile in bamlist:
	barcode = 'IonXpress_' + bamfile.split('IonXpress_')[-1].split('.bam')[0]
	barcode2bam[barcode] = bamfile
	sample = bamfile.split('/')[-1].split('_IonXpress')[0]
	barcode2sample[barcode] = sample
	target = barcodes_json[barcode]['target_region_filepath'].split('/')[-1]
	for _run_type in global_param['run_type']:
		if global_param['run_type'][_run_type]['target_bed'].split('/')[-1] == target:
			barcode2runtype[barcode] = _run_type
			break

if not (('LAM' in barcode2runtype.values()) or ('LAM_2018' in barcode2runtype.values())):
	exit() # no LAM patient, no dupG
if not os.path.isdir(dupg_dir):
	subprocess.call(['mkdir',dupg_dir])

#sample2data = {}
for barcode in barcode2runtype:
	sample = barcode2sample[barcode]
	sample = sample.replace(' ','_')
	iscontrol = False # check if barcode is not a control
	for control_name in control_names:
		if control_name in sample.upper():
			iscontrol = True
	if iscontrol:
		continue
	if barcode2runtype[barcode] not in ['LAM','LAM_2018']:
		continue
		
	bamfile = pysam.AlignmentFile(barcode2bam[barcode],'rb')
	filtered_bam = pysam.AlignmentFile('%s/dupG/HotCount/fastq_temp/%s_%s.bam' % (pipeline_folder,sample,barcode), 'wb', template=bamfile)
	for read in bamfile.fetch('chr20',31022430,31022460):
		filtered_bam.write(read)

###############################
#######  get checkMut   #######
###############################

# checkMut already done via routine
try:
	subprocess.call(['cp', options.run_folder+'/_checkMut/ASXL1_c.1927dup.png', dupg_dir+'/ASXL1_c.1927dup.png'])
	subprocess.call(['cp', options.run_folder+'/_checkMut/ASXL1_c.1927del.png', dupg_dir+'/ASXL1_c.1927del.png'])
except:
	pass

###############################
####### HotCount script #######
###############################

print "HotCount script"
# generer les fastq necessaires pour le script
for barcode in barcode2runtype:
	sample = barcode2sample[barcode]
	sample = sample.replace(' ','_')
	iscontrol = False # check if barcode is not a control
	for control_name in control_names:
		if control_name in sample.upper():
			iscontrol = True
	if iscontrol:
		continue
	if barcode2runtype[barcode] not in ['LAM','LAM_2018']:
		continue
	subprocess.call(['samtools','fastq','%s/dupG/HotCount/fastq_temp/%s_%s.bam' % (pipeline_folder,sample,barcode)],stdout=open('%s/dupG/HotCount/fastq_temp/%s_%s.fastq' % (pipeline_folder,sample.replace(' ','@'),barcode),'w'))

# lancement script
subprocess.call(['bash','%s/dupG/HotCount/do_it.sh' % pipeline_folder,'%s/dupG/HotCount/asxl1.txt' % pipeline_folder,'%s/dupG/HotCount/fastq_temp/*.fastq' % pipeline_folder],stdout=open('%s/dupG/HotCount/fastq_temp/comtpage_asxl1.tsv' % pipeline_folder,'w'))
subprocess.call(['Rscript','%s/dupG/HotCount/do_it.R' % pipeline_folder,'%s/dupG/HotCount/fastq_temp/comtpage_asxl1.tsv' % pipeline_folder,'ALL','dupG','delG'],stdout=open('%s/dupG/HotCount/fastq_temp/HotCount_results.tsv' % pipeline_folder,'w'))

# copier resultat
subprocess.call(['cp','%s/dupG/HotCount/fastq_temp/HotCount_results.tsv' % pipeline_folder,dupg_dir])

# nettoyer 
subprocess.call(['rm %s/dupG/HotCount/fastq_temp/*.*' % pipeline_folder],shell=True)
