#!/usr/bin/python
import subprocess
import glob
import sys
import time
import datetime
import os

def update_completed_runs_list(completed_run_list_path,completed_runs):
	with open(completed_run_list_path,'w') as crl:
		crl.writelines(completed_runs)

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
ngs_folder = '/media/n06lbth/sauvegardes_pgm'
illumina_folder_path = '/media/n06lbth/sauvegardes_pgm/Illumina_output'
completed_run_list_path = '%s/completed_run_list.txt' % illumina_folder_path

illumina_folders = glob.glob('%s/*/' % illumina_folder_path)
illumina_folders = [folder.split('/')[-2] for folder in illumina_folders]

completed_runs = {}
with open(completed_run_list_path,'r') as crl:
	for line in crl:
		line = line.replace('\n','').replace('\r','')
		line = line.split('\t')
		completed_runs[line[0]] = line[1]

new_analysis = False
# IF NEW RUN
for folder in illumina_folders:
	if folder not in completed_runs.keys():
		new_analysis = True
		new_folder = '%s/%s' % (illumina_folder_path,folder)
		print "NEW FOLDER DETECTED : %s" % new_folder
		copy_complete_file = '%s/CopyComplete.txt' % new_folder

		completed_runs[folder] = 'waiting CopyComplete.txt...'
		update_completed_runs_list(completed_run_list_path,completed_runs)
		while not os.path.isfile(copy_complete_file):
			print "waiting for CopyComplete.txt..."
			time.sleep(600) # try every 10min

		print "- launching prepare_illumina_run_folder.py ..."
		completed_runs[folder] = 'bcl2fastq...'
		update_completed_runs_list(completed_run_list_path,completed_runs)
		subprocess.call(['python','%s/prepare_illumina_run_folder.py' % pipeline_folder,'--illumina_folder',new_folder])

		############################################################################################################################

		print "- parsing SampleSheet ..."
		run_project = ''
		sub_project = ''
		sample_sheet = '%s/SampleSheet.csv' % new_folder
		ss_reader = open(sample_sheet,'r')
		experiment_name = 'noname_run_%s' % time.strftime("%d-%m-%Y-%Hh%M")

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
		completed_runs[folder] = 'run_analysis in progress...'
		update_completed_runs_list(completed_run_list_path,completed_runs)
		subprocess.call(['python','%s/run_analysis.cross.py' % pipeline_folder,'--run',run_folder])
		
		print "- run analysis is completed"
		completed_runs[folder] = 'completed'
		update_completed_runs_list(completed_run_list_path,completed_runs)

if not new_analysis:
	print "[%s] nothing to do" % datetime.date.today()
