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
import urllib2

############################################################################################
parser = OptionParser()
parser.add_option('-r', '--run-folder', help='Run folder ', 								dest='run_folder')
parser.add_option('-c', '--chrom', 		help='chromosome', 									dest='chrom')
parser.add_option('-s', '--start', 		help='start pos to search', 						dest='start')
parser.add_option('-e', '--end', 		help='end pos to search', 							dest='end')
parser.add_option('-w', '--wt', 		help='wild-type sequence', 							dest='wt')
parser.add_option('-v', '--variant', 	help='variant sequence', 							dest='variant')
parser.add_option('-n', '--mutname', 	help='short description of the searched variant', 	dest='mutname')
parser.add_option('-y', '--nm', 		help='Alternative method with NM (and c_pos)', 		dest='nm', default=False)
parser.add_option('-z', '--cpos', 		help='Alternative method with c_pos (and NM)', 		dest='cpos', default=False)
parser.add_option('-x', '--reverse', 	help='reverse sequence from mutalyzer', 			dest='reverse', default=False, action='store_true')
parser.add_option('-t', '--run-type', 	help='filter sample by run type', 					dest='runtype', default=False)
parser.add_option('-a', '--wsize', 		help='larger windows size', 						dest='wsize', default=3)
(options, args) = parser.parse_args()

wsize = int(options.wsize)
with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)
with open(global_param['NM_data'], 'r') as nmdata:
	NM_version = json.load(nmdata)

if os.path.isfile(options.run_folder+'/barcodes.json'):
	with open(options.run_folder+'/barcodes.json', 'r') as g:
		barcodes_json = json.load(g)
else:
	print "error : barcodes.json not found in run folder"

bamlist = glob.glob(options.run_folder+'/*/*.bam')
bamlist = [item for item in bamlist if not 'processed' in item]

barcode2sample = {}
barcode2runtype = {}
barcode2bam = {}
for bamfile in bamlist:
	barcode = 'IonXpress_' + bamfile.split('IonXpress_')[-1].split('.bam')[0]
	barcode2bam[barcode] = bamfile
	sample = bamfile.split('/')[-1].split('_IonXpress')[0]
	iscontrol = False
	for control_name in ['NTC_','H2O']:
		if control_name in sample.upper():
			iscontrol = True
			break
	if iscontrol:
		continue
	barcode2sample[barcode] = sample
	target = barcodes_json[barcode]['target_region_filepath'].split('/')[-1]
	for _run_type in global_param['run_type']:
		if global_param['run_type'][_run_type]['target_bed'].split('/')[-1] == target:
			barcode2runtype[barcode] = _run_type
			break

if options.runtype:
	if options.runtype not in barcode2runtype.values():
		print "- CheckMut : runtype not found, exit"
		sys.exit()
		
checkMut = options.run_folder + '/_checkMut'
if not os.path.isdir(checkMut):
	subprocess.call(['mkdir',checkMut])

if options.nm:
	nm = options.nm + '.' + NM_version[options.nm]['version']
	cpos = options.cpos
	cpos = cpos.replace('+','%2B')  # e.g : c.1959+1G>A
	cpos = cpos.replace('-','%2D')  # e.g : c.1912-2A>C
	url_mutalyzer = "https://mutalyzer.nl/json/runMutalyzer?variant=%s:%s" % (nm,cpos)
	f = urllib2.urlopen(url_mutalyzer)
	d = json.loads(f.read())
	genename = d['legend'][0]['name']
	genename = genename.split('_')[0]
	print "Gene : " + genename
	visu = d['rawVariants'][0]['visualisation']
	
	visu1 = visu.split('\n')[0]
	leftseq1 = visu1.split(' ')[0][-wsize:]
	middleseq1 = visu1.split(' ')[1]
	rightseq1 = visu1.split(' ')[2][:wsize]
	
	visu2 = visu.split('\n')[1]
	leftseq2 = visu2.split(' ')[0][-wsize:]
	middleseq2 = visu2.split(' ')[1]
	rightseq2 = visu2.split(' ')[2][:wsize]
	
	WT = leftseq1 + middleseq1 + rightseq1
	VARIANT = leftseq2 + middleseq2 + rightseq2
	
	if options.reverse:
		WT = WT.replace('A','W').replace('T','X').replace('G','Y').replace('C','Z')
		VARIANT = VARIANT.replace('A','W').replace('T','X').replace('G','Y').replace('C','Z')
		WT = WT.replace('W','T').replace('X','A').replace('Y','C').replace('Z','G')
		VARIANT = VARIANT.replace('W','T').replace('X','A').replace('Y','C').replace('Z','G')
		WT = WT[::-1]
		VARIANT = VARIANT[::-1]
	
	WT = WT.replace('-','')
	VARIANT = VARIANT.replace('-','')
	
	for v in range(int(NM_version[nm]['version']),0,-1):
		try:
			nmv = nm + '.' + str(int(v))
			url_mutalyzer = "https://mutalyzer.nl/json/numberConversion?build=hg19;variant=%s:%s" % (nmv,cpos)
			f = urllib2.urlopen(url_mutalyzer)
			d = json.loads(f.read())
			chrom = d[0].split('NC_')[-1].split('.')[0]
		except:
			continue
		break
	if (int(chrom) == 23):
		chrom = 'chrX'
	elif (int(chrom) == 24):
		chrom = 'chrY'
	else:
		chrom = 'chr' + str(int(chrom))
	genopos = d[0].split('g.')[-1].split('_')[0].split('A')[0].split('T')[0].split('G')[0].split('C')[0].split('del')[0].split('dup')[0]
	start_search = int(genopos)-1
	stop_search = int(genopos)+1
	mutname = genename + '_' + cpos.replace('>','_').replace('<','_')
	print "%s:%s-%s" % (chrom,start_search,stop_search)
	print "%s  ->  %s" % (WT,VARIANT)
	print "results for %s" % mutname
	
else:
	chrom = options.chrom 			# 'chr17'
	start_search = int(options.start) 	# 37880973
	stop_search = int(options.end) 		# 37881000
	WT = options.wt 				# "GAAGCATACGTGATGGCT"
	VARIANT = options.variant 		# "GAAGCATACGTGATGGCATACGTGATGGCT"
	mutname = options.mutname 		# 'ERBB2_2313_2324dup_COSM20959'

############################################################################################

sample2data = {}
for barcode in barcode2runtype:
	if options.runtype:
		if barcode2runtype[barcode] != options.runtype:
			continue
	sample = barcode2sample[barcode]
	sample = sample.replace(' ','_')
	depth = 0
	WT_count = 0
	variant_count = 0
	try:
		bamfile = pysam.AlignmentFile(barcode2bam[barcode],'rb')
		for read in bamfile.fetch(chrom,start_search,stop_search):
			depth = depth +1
			seq = read.query_sequence
			if WT in seq :
				WT_count = WT_count+1
			if VARIANT in seq :
				variant_count = variant_count+1
	except:
		continue

	sample2data[sample] = {}
	sample2data[sample]['DP'] = depth
	sample2data[sample]['VARIANT'] = variant_count
	sample2data[sample]['WT'] = WT_count

### CALCUL DE LA DETECTION DE DUPG ###
list_freq = []
for sample in sample2data:
	if sample2data[sample]['DP'] == 0:
		freq = 0.0
	else:
		freq = numpy.divide(float(sample2data[sample]['VARIANT']),float(sample2data[sample]['DP']))*100
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

csvfile = open('%s/%s.csv' % (checkMut,mutname),'wb')
results = csv.writer(csvfile, delimiter =';')
#header
results.writerow([mutname,chrom,'search from ' + str(start_search),'to ' + str(stop_search),'wt seq : ' + WT,'variant seq : ' + VARIANT])
results.writerow([mutname,'mean = %s' % numpy.around(moy_freq,1),'standard_deviation = %s' % numpy.around(ecart_type_freq,1), 'threshold = %s' % numpy.around(threshold,1)])
results.writerow(['Sample',mutname,'Depth','Freq VARIANT','is_real'])
for sample in sample2data.keys():	
	results.writerow([sample,sample2data[sample]['VARIANT'],sample2data[sample]['DP'],numpy.around(sample2data[sample]['FREQ'],1),sample2data[sample]['ISREAL']])
	print "- {0:40} {1}".format(sample,str(numpy.around(sample2data[sample]['FREQ'],1)))
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
ax.set_ylabel('VARIANT FREQ',fontsize=12)
xticks(ind,sample2data.keys(),fontsize=4)
ax.plot([0., N], [threshold, threshold], 'k--', color='blue') # horizontal line indicating the threshold
max_y = max(tab_freq) #ax.set_xlim([0,max(tab_freq_del)])
ax.set_ylim([0,max_y*1.5])	#augmentation un peu de l'ordonne, trop courte
fig.autofmt_xdate()
fig.savefig('%s/%s.png' % (checkMut,mutname),bbox_inches='tight',dpi=600)

