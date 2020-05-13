#!/usr/bin/env python
import os
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

"""Script searching EGFR p.(Ile789Phe) mutation"""

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

mut_dir = options.run_folder + '/_I789F_search'
WT = "GCTCATCACGCA"
MUT = "GCTCTTCACGCA"
mutname = 'I789F'
chromosome = 'chr7'
start_pos_search = 55249052 # pas besoin d'etre precis, pour faire une fenetre de recherche
end_pos_search = 55249092

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

#if not 'SBT' in barcode2runtype.values():
	#exit()
if not os.path.isdir(mut_dir):
	subprocess.call(['mkdir',mut_dir])

sample2data = {}
for barcode in barcode2runtype:
	sample = barcode2sample[barcode]
	sample = sample.replace(' ','_')
	depth = 0
	WT_count = 0
	mut_count = 0
	sample2data[sample] = {}
	bamfile = pysam.AlignmentFile(barcode2bam[barcode],'rb')
	for read in bamfile.fetch(chromosome,start_pos_search,end_pos_search):
		depth = depth +1
		seq = read.query_sequence
		if WT in seq :
			WT_count = WT_count+1
		if MUT in seq :
			mut_count = mut_count+1
	
	sample2data[sample]['DP'] = depth
	sample2data[sample]['MUT'] = mut_count
	sample2data[sample]['WT'] = WT_count

### CALCUL DE LA DETECTION DE DUPG ###
list_freq = []
for sample in sample2data:
	if sample2data[sample]['DP'] == 0:
		freq = 0.0
	else:
		freq = numpy.divide(float(sample2data[sample]['MUT']),float(sample2data[sample]['DP']))*100
	sample2data[sample]['FREQ'] = freq
	list_freq.append(freq)
	
moy_freq = numpy.mean(list_freq)
ecart_type_freq = numpy.std(list_freq,ddof=1)
threshold = moy_freq + 2*ecart_type_freq

for sample in sample2data.keys():
	if sample2data[sample]['FREQ'] > threshold:
		sample2data[sample]['ISREAL'] = "Probable"
	else:
		sample2data[sample]['ISREAL'] = "Non"

csvfile = open('%s/%s.csv' % (mut_dir,mutname),'w')
results = csv.writer(csvfile, delimiter =';')
results.writerow([mutname,'mean = %s' % numpy.around(moy_freq,1),'standard_deviation = %s' % numpy.around(ecart_type_freq,1), 'threshold = %s' % numpy.around(threshold,1)])
results.writerow(['Sample',mutname,'Depth','Freq MUT','is_real'])
for sample in sample2data.keys():	
	results.writerow([options.run_folder.split('/')[-1],sample,sample2data[sample]['MUT'],sample2data[sample]['DP'],numpy.around(sample2data[sample]['FREQ'],1),sample2data[sample]['ISREAL']])
csvfile.close()	
