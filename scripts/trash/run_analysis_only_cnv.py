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
	
	bam_data[bamfile] = {'barcode':barcode,'sample':sample,'run_type':run_type,'reference':reference,'sample_folder':sample_folder,'intermediate_folder':intermediate_folder,'target_bed':target_bed,'bed_merged':bed_merged,'param':param,'param_hotspot_only':param_hotspot_only,'hotspot_vcf':hotspot_vcf}
		
#######################
## COVERAGE ANALYSIS ## 
#######################

print " [%s] Coverage Analysis ..." % (time.strftime("%H:%M:%S"))	
for bamfile in bam_data:
	coverage_folder = '%s/coverage' % bam_data[bamfile]['intermediate_folder']
	if not os.path.isdir(coverage_folder):
		subprocess.call(['mkdir', coverage_folder])
		
	cmd = subprocess.Popen([
		'bash', '/DATA/work/coverageAnalysis/run_coverage_analysis.sh',
		'-ag',
		'-D',coverage_folder,
		'-B',bam_data[bamfile]['target_bed'],
		bam_data[bamfile]['reference'],
		bamfile
		], stdout=open(coverage_folder+'/run_coverage_analysis.stdout.txt','w'), stderr=open(coverage_folder+'/run_coverage_analysis.stderr.txt','w'))	
		#		'-L', 'hg19',
	cmd.communicate()
	
time.sleep(60)

##########
### CNV ##
##########

if options.full_run:
	if not os.path.isdir(run_folder+'/_CNA'):
		subprocess.call(['mkdir', run_folder+'/_CNA'])
	print " [%s] CNV Analysis ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen(['python','/DATA/work/CNV/run_cna.py', run_folder], stdout=open(run_folder+'/_CNA/run_cna.stdout.txt','w'), stderr=open(run_folder+'/_CNA/run_cna.stderr.txt','w'))
	cmd.communicate()
####

time.sleep(60)

for bamfile in bam_data:
	shutil.rmtree(bam_data[bamfile]['intermediate_folder'])
