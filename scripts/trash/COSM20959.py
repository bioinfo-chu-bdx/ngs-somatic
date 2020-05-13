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

"""Script searching COSM20959 mutation"""

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

#control_names = ['H2O','H20','NTC'] # liste des noms possibles pour les temoins negatifs
checkMut = options.run_folder + '/_checkMut'
WT = "GAAGCATACGTGATGGCT"
COSM20959 = "GAAGCATACGTGATGGCATACGTGATGGCT"

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
if not os.path.isdir(checkMut):
	subprocess.call(['mkdir',checkMut])

sample2data = {}
for barcode in barcode2runtype:
	sample = barcode2sample[barcode]
	sample = sample.replace(' ','_')
	#iscontrol = False # check if barcode is not a control
	#for control_name in control_names:
		#if control_name in sample.upper():
			#iscontrol = True
	#if iscontrol:
		#continue
	#if barcode2runtype[barcode] != 'SBT':
		#continue
		
	#print "Processing sample : " + sample
	depth = 0
	WT_count = 0
	cosm20959_count = 0
	
	sample2data[sample] = {}
	
	bamfile = pysam.AlignmentFile(barcode2bam[barcode],'rb')
	for read in bamfile.fetch('chr17',37880973,37881000):
		depth = depth +1
		seq = read.query_sequence
		if WT in seq :
			WT_count = WT_count+1
		if COSM20959 in seq :
			cosm20959_count = cosm20959_count+1
	
	sample2data[sample]['DP'] = depth
	sample2data[sample]['COSM20959'] = cosm20959_count
	sample2data[sample]['WT'] = WT_count

### CALCUL DE LA DETECTION DE DUPG ###
list_freq = []
for sample in sample2data:
	if sample2data[sample]['DP'] == 0:
		freq = 0.0
	else:
		freq = numpy.divide(float(sample2data[sample]['COSM20959']),float(sample2data[sample]['DP']))*100
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

csvfile = open('%s/ERBB2_2313_2324dup_COSM20959.csv' % checkMut,'wb')
results = csv.writer(csvfile, delimiter =';')
#header
results.writerow(['cosm20959','mean = %s' % numpy.around(moy_freq,1),'standard_deviation = %s' % numpy.around(ecart_type_freq,1), 'threshold = %s' % numpy.around(threshold,1)])
results.writerow(['Sample','cosm20959','Depth','Freq COSM20959','is_real'])
for sample in sample2data.keys():	
	results.writerow([sample,sample2data[sample]['COSM20959'],sample2data[sample]['DP'],numpy.around(sample2data[sample]['FREQ'],1),sample2data[sample]['ISREAL']])
csvfile.close()	

### GRAPH ###
fig = figure()

N = len(sample2data.keys())
ind = numpy.arange(N)
width = 0.4

tab_freq = []
for sample in sample2data.keys():
	tab_freq.append(sample2data[sample]['FREQ'])

ax = fig.add_subplot(111)
rects = ax.bar(ind, tab_freq, width, facecolor='#9999ff', label='d')
ax.set_xlabel('Patients',fontsize=12)
ax.set_ylabel('COSM20959 FREQ',fontsize=12)
xticks(ind,sample2data.keys(),fontsize=4)
ax.plot([0., N], [threshold, threshold], 'k--', color='blue') # horizontal line indicating the threshold
max_y = max(tab_freq) #ax.set_xlim([0,max(tab_freq_del)])
ax.set_ylim([0,max_y*1.5])	#augmentation un peu de l'ordonne, trop courte
fig.autofmt_xdate()
fig.savefig('%s/ERBB2_2313_2324dup_COSM20959.png' % checkMut,bbox_inches='tight',dpi=600)

