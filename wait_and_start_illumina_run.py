#!/usr/bin/python
import subprocess
import sys
import time
import os

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
ngs_folder = '/media/n06lbth/sauvegardes_pgm'

illumina_output_folder = sys.argv[1]
copy_complete_file = '%s/CopyComplete.txt' % illumina_output_folder

while not os.path.isfile(copy_complete_file):
	print "waiting for CopyComplete.txt..."
	time.sleep(600) # try every 10min

print "- launching prepare_illumina_run_folder.py ..."
subprocess.call(['python','%s/prepare_illumina_run_folder.py' % pipeline_folder,'--illumina_folder',illumina_output_folder])

############################################################################################################################

print "- parsing SampleSheet ..."
run_project = ''
sub_project = ''
sample_sheet = '%s/SampleSheet.csv' % illumina_output_folder
ss_reader = open(sample_sheet,'r')
experiment_name = 'unnamed_run_%s' % time.strftime("%d-%m-%Y-%Hh%M")

for line in ss_reader:
	if line in ['\n','\r\n']:
		continue
	if line.startswith('Experiment Name'):
		experiment_name = line.replace('\n','').replace('\r','').split(',')[-1]
	if line.startswith('Run Project'):
		run_project = line.replace('\n','').replace('\r','').split(',')[-1]
	if line.startswith('Sub Project'):
		sub_project = line.replace('\n','').replace('\r','').split(',')[-1]

if sub_project == '':
	output_location = '%s/%s' % (ngs_folder,run_project)
else:
	output_location = '%s/%s/%s' % (ngs_folder,run_project,sub_project)
run_folder = '%s/%s' % (output_location,experiment_name)

print "- launching run_analysis.cross.py ..."
subprocess.call(['python','%s/run_analysis.cross.py' % pipeline_folder,'--run',run_folder])

