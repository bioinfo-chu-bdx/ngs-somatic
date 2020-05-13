#!/usr/bin/python
import sys
import os
import time
import glob
import json
import subprocess
import zipfile
import shutil
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
fp_file = global_param['fp_file']
control_names = ['H2O','H20','NTC'] # liste des noms possibles pour les temoins negatifs
checkconta_bamlist = []
sampleask = False

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
	if zipfile.is_zipfile(intermediate_folder+'.zip'):
		z = zipfile.ZipFile(intermediate_folder+'.zip','r')
		z.extractall(intermediate_folder)
		z.close()
	
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

	#################
	## FINALREPORT ##
	#################
	
	if not os.path.isdir(intermediate_folder+'/finalReport'):
		subprocess.call(['mkdir', intermediate_folder+'/finalReport'])
	finalreport_name = sample+'_'+barcode+'.finalReport.xlsx'
	finalReport_path = intermediate_folder+'/finalReport/'+finalreport_name
	
	print "\t - [%s] finalReport.py ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen([
		'python','/DATA/work/finalReport/finalReport.py',
		'--bam', bamfile,
		'--run-type', run_type,
		'--output', finalReport_path,
	], stdout=open(intermediate_folder+'/finalReport/finalReport.stdout.txt','w'), stderr=open(intermediate_folder+'/finalReport/finalReport.stderr.txt','w'))
	cmd.communicate()

	################################################

	subprocess.call(['mv',finalReport_path,sample_folder])
	subprocess.call(['mv',intermediate_folder+'/_ALAMUT_load_variants.vbs',sample_folder])
	subprocess.call(['mv',intermediate_folder+'/_ALAMUT_load_processed_bam.vbs',sample_folder])
	subprocess.call(['mv',intermediate_folder+'/_PRINT_finalReport.vbs',sample_folder])
	
	shutil.make_archive(intermediate_folder,'zip',intermediate_folder)
	shutil.rmtree(intermediate_folder)

	print " [%s] ... Done" % (time.strftime("%H:%M:%S"))
	
	
if options.full_run and os.path.isfile(run_folder+'/barcodes.json'):
	
	### VBS scripts for printing
	print " [%s] Making VBS for fast printing..." % (time.strftime("%H:%M:%S"))
	subprocess.call(['python','/DATA/work/finalReport/make_vbs_print.py',run_folder,run_folder+'/barcodes.json'])
	
