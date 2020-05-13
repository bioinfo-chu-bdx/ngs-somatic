#!/usr/bin/python
import sys
import os
import glob
import subprocess
import json
import csv
import zipfile
import shutil
import time
import numpy
from optparse import OptionParser

### GATHERING PARAMETERS ############################################################

parser = OptionParser()
parser.add_option('-b', '--bam',		help="Input bam file for SINGLE BAM ANALYSIS (must contain flow signal from TSS)", dest='bam')
parser.add_option('-r', '--full-run',	help="Run folder path  for FULL RUN ANALYSIS", dest='full_run') 
parser.add_option('-t', '--run-type',	help="Run type : SBT, LAM, TP53-HEMATO... (OPTIONAL, NOT NEEDED if barcodes.json is present in run folder)", dest='run_type') 
(options, args) = parser.parse_args()

if options.bam and options.full_run:
	sys.stderr.write("[run_analysis.py] Error: <--bam> and <--full-run> are not compatibles\n")
	sys.exit()
if options.full_run:
	bamlist = glob.glob(options.full_run+'/*/*.bam')
	bamlist = [item for item in bamlist if not 'processed' in item]
elif options.bam:
	bamlist = [options.bam]
else:
	sys.stderr.write("[run_analysis.py] Error: no <--bam> or <--full-run> specified\n")
	sys.exit()
	
with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)
sampleask = True

### RUN ANALYSIS ######################################################################

bam_data = {}
	
for bamfile in sorted(bamlist) :
	sample = bamfile.split('/')[-1].split('_IonXpress')[0]
	barcode = 'IonXpress_' + bamfile.split('IonXpress_')[-1].split('.bam')[0]
	sample_folder = os.path.dirname(bamfile)
	run_folder = os.path.dirname(sample_folder)
	
	### RUN TYPE ###	
	run_type = None
	# if run-type manually specified
	if options.run_type:
		run_type = options.run_type
		if run_type not in global_param['run_type'].keys():
			sys.stderr.write("[run_analysis.py] Error: run type <%s> not found in global_parameters.json file\n" % run_type)
			sys.exit()
	# default : via barcodes.json in run folder
	elif os.path.isfile(run_folder+'/barcodes.json'):
		with open(run_folder+'/barcodes.json', 'r') as g:
			barcodes_json = json.load(g)
			bed_name = barcodes_json[barcode]['target_region_filepath'].split('/unmerged/detail/')[-1]
			for _run_type in global_param['run_type']:
				if global_param['run_type'][_run_type]['target_bed'].split('/')[-1] == bed_name:
					run_type = _run_type
					break
					
	### RUN TYPE PARAMETERS ###	
	if not run_type:
		print "\t -- Error : run type not found. Sample %s will not be processed." % sample
		continue
	reference = global_param['run_type'][run_type]['reference']
	target_bed = global_param['run_type'][run_type]['target_bed']
	bed_merged = global_param['run_type'][run_type]['merged_bed']
	param = global_param['run_type'][run_type]['vc_parameters']
	param_hotspot_only = global_param['run_type'][run_type]['vc_parameters_hotspot_only']
	hotspot_vcf = global_param['run_type'][run_type]['hotspot_vcf']

	# dossier results
	intermediate_folder = sample_folder + '/intermediate_files'
	if not os.path.isdir(intermediate_folder):
		subprocess.call(['mkdir', intermediate_folder])
	#if zipfile.is_zipfile(intermediate_folder+'.zip'):
		#z = zipfile.ZipFile(intermediate_folder+'.zip','r')
		#z.extractall(intermediate_folder)
		#z.close()
	
	bam_data[bamfile] = {'barcode':barcode,'sample':sample,'run_type':run_type,'reference':reference,'sample_folder':sample_folder,'intermediate_folder':intermediate_folder,'target_bed':target_bed,'bed_merged':bed_merged,'param':param,'param_hotspot_only':param_hotspot_only,'hotspot_vcf':hotspot_vcf}

### RESTE DE L'ANALYSE PATIENT PAR PATIENT ###	
for bamfile in sorted(bamlist) :
	print " [%s] Processing %s ..." % (time.strftime("%H:%M:%S"),bamfile.split('/')[-1])
	
	###################################################################
	# FOR SKIPING SAMPLES # WARNING check-contamination will not work #
	if sampleask:
		proceed = raw_input('continue? (y/n/stopask)\n')              
		if proceed == 'y' or proceed == 'yes':                        
			pass
		elif proceed == 'stopask':
			sampleask = False
			pass                                                       
		else:                                                         
			continue                                                  
	###################################################################
	
	print "\t * panel %s" % run_type	
	barcode = bam_data[bamfile]['barcode']
	sample = bam_data[bamfile]['sample']
	run_type = bam_data[bamfile]['run_type']
	reference = bam_data[bamfile]['reference']
	sample_folder = bam_data[bamfile]['sample_folder']
	intermediate_folder = bam_data[bamfile]['intermediate_folder']
	target_bed = bam_data[bamfile]['target_bed']
	bed_merged = bam_data[bamfile]['bed_merged']
	param = bam_data[bamfile]['param']
	param_hotspot_only = bam_data[bamfile]['param_hotspot_only']
	hotspot_vcf = bam_data[bamfile]['hotspot_vcf']

	######################
	## COUNT BASES 100X ##
	######################
	
	depth_file = open(intermediate_folder+'/tvc_de_novo/depth.txt','r')
	depth_reader = csv.reader(depth_file,delimiter='\t')

	total = 0
	sup100X = 0
	depths = []
	for line in depth_reader:
		total = total + 1
		depth = int(line[2])
		depths.append(depth)
		if depth >= 100:
			sup100X = sup100X +1
		
	meandepth = numpy.mean(depths)
	print "--> %s percent bases >100X" % str(float(sup100X)/float(total))
	print "--> mean reads per base = %s" % str(meandepth)
	
	#shutil.make_archive(intermediate_folder,'zip',intermediate_folder)
	#shutil.rmtree(intermediate_folder)
