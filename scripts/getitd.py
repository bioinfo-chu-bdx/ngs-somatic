#!/usr/bin/python
import os
import sys
import time
import glob
import json
import subprocess
from datetime import date
from optparse import OptionParser

def normalize_folder_path(folder_path): # REMOVE TRAILING SLASH /
	if os.path.isdir(folder_path):
		if folder_path.endswith('/'):
			folder_path = folder_path[:-1]
	return folder_path

### GATHERING PARAMETERS ############################################################

FNULL = open(os.devnull, 'w')
parser = OptionParser()
parser.add_option('-r', '--run',					help="run folder path for FULL RUN ANALYSIS",dest='run') 
parser.add_option('-s', '--sample',					help="sample folder path for SINGLE SAMPLE ANALYSIS", dest='sample')
(options, args) = parser.parse_args()

if options.run and options.sample:
	sys.stderr.write("[run_analysis.py] Error: choose either <--run> or <--sample>, not both\n")
	sys.exit()
if options.run:
	options.run = normalize_folder_path(options.run)
	run_folder = options.run
elif options.sample:
	options.sample = normalize_folder_path(options.sample)
	sample_folder = options.sample
	run_folder = os.path.dirname(sample_folder)
else:
	sys.stderr.write("[run_analysis.py] Error: no <--fastq> or <--full-run> specified\n")
	sys.exit()

if os.path.isfile('%s/barcodes.json' % run_folder):
	with open('%s/barcodes.json' % run_folder, 'r') as g:
		barcodes_json = json.load(g)
else:
	sys.exit("barcodes.json file not found")

if options.sample:
	for barcode in barcodes_json.keys():
		if barcodes_json[barcode]['sample'] != os.path.basename(options.sample):
			del barcodes_json[barcode]
else:
	for barcode in barcodes_json.keys():
		if 'checkContamination' in barcodes_json[barcode]['sample']:
			del barcodes_json[barcode]

ordered_barcodes = [item[1] for item in sorted([(barcodes_json[barcode]['sample'],barcode) for barcode in barcodes_json])]

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))


getitd_folder = '%s/_getitd' % run_folder # pour les stdout et stderr
if not os.path.isdir(getitd_folder):
	subprocess.call(['mkdir', getitd_folder])

ps_list = []
processing = True
os.chdir('%s/getITD/' % pipeline_folder)
for barcode in ordered_barcodes:
	barcodes_json[barcode]['bam'] = '%s/%s/%s_%s.bam' % (run_folder,barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode)
	barcodes_json[barcode]['flt3_bam'] = '%s/%s_%s.flt3-filtered.bam' % (getitd_folder,barcodes_json[barcode]['sample'],barcode)
	barcodes_json[barcode]['flt3_fastq'] = '%s/%s_%s.flt3-filtered.fastq' % (getitd_folder,barcodes_json[barcode]['sample'],barcode)
	print " [%s] - %s : bedtools intersect BAM - FLT3 ..." % (barcodes_json[barcode]['sample'],time.strftime("%H:%M:%S"))
	ps = subprocess.Popen(['bedtools','intersect','-a',barcodes_json[barcode]['bam'],'-b','%s/reference_files/FLT3.bed' % pipeline_folder],stdout=open(barcodes_json[barcode]['flt3_bam'],'w'))
	ps_list.append(ps)

while processing :
	ps_states = [ps.poll() is not None for ps in ps_list]
	if False not in ps_states :
		processing = False

ps_list = []
processing = True
for barcode in ordered_barcodes:
	print " [%s] - %s : bedtools bamtofastq ..." % (barcodes_json[barcode]['sample'],time.strftime("%H:%M:%S"))
	ps = subprocess.Popen(['bedtools','bamtofastq','-i',barcodes_json[barcode]['flt3_bam'],'-fq',barcodes_json[barcode]['flt3_fastq']])
	ps_list.append(ps)

while processing :
	ps_states = [ps.poll() is not None for ps in ps_list]
	if False not in ps_states :
		processing = False

ps_list = []
processing = True
for barcode in ordered_barcodes:
	print " [%s] - %s : getitd.py ..." % (barcodes_json[barcode]['sample'],time.strftime("%H:%M:%S"))
	ps = subprocess.Popen(['python3','getitd.py','%s_%s' % (barcodes_json[barcode]['sample'],barcode),barcodes_json[barcode]['flt3_fastq']])
	ps_list.append(ps)

while processing :
	ps_states = [ps.poll() is not None for ps in ps_list]
	if False not in ps_states :
		processing = False

ps_list = []
processing = True
for barcode in ordered_barcodes:
	print " [%s] - %s : move / delete temp files ..." % (barcodes_json[barcode]['sample'],time.strftime("%H:%M:%S"))
	subprocess.call(['mv', '%s_%s_getitd' % (barcodes_json[barcode]['sample'],barcode), getitd_folder])
	subprocess.call(['rm', barcodes_json[barcode]['flt3_bam']])
	subprocess.call(['rm', barcodes_json[barcode]['flt3_fastq']])
