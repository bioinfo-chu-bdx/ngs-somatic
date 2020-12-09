#!/usr/bin/python
import subprocess
import glob
import time
import os

def update_completed_runs_list(completed_run_list_path,completed_runs):
	with open(completed_run_list_path,'w') as crl:
		for completed_run in sorted(completed_runs.keys()):
			line = '%s\t%s\n' % (completed_run,completed_runs[completed_run])
			crl.write(line)

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
		print "[%s] New folder detected : %s" % (time.strftime("%Y/%m/%d]  [%H:%M"),new_folder)
		copy_complete_file = '%s/CopyComplete.txt' % new_folder

		# CHECK IF FOLDER IS AN ILLUMINA OUTPUT FOLDER
		if not os.path.isfile('%s/RunParameters.xml' % new_folder):
			print "[%s] folder does not appear to be an illumina output folder, skipping" % time.strftime("%Y/%m/%d]  [%H:%M")
			completed_runs[folder] = 'not illumina'
			update_completed_runs_list(completed_run_list_path,completed_runs)

		completed_runs[folder] = 'waiting CopyComplete.txt...'
		update_completed_runs_list(completed_run_list_path,completed_runs)
		while not os.path.isfile(copy_complete_file):
			print "[%s] waiting for CopyComplete.txt..." % time.strftime("%Y/%m/%d]  [%H:%M")
			time.sleep(600) # try every 10min

		print "[%s] - launching prepare_illumina_run_folder.py ..." % time.strftime("%Y/%m/%d]  [%H:%M")
		completed_runs[folder] = 'bcl2fastq...'
		update_completed_runs_list(completed_run_list_path,completed_runs)
		cmd = subprocess.Popen(['python','%s/prepare_illumina_run_folder.py' % pipeline_folder,'--illumina_folder',new_folder],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		cmd.communicate()

		############################################################################################################################

		print "[%s] - parsing SampleSheet ..." % time.strftime("%Y/%m/%d]  [%H:%M")
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

		if run_project == 'FusionPlex_CTL':
			completed_runs[folder] = 'completed (only FastQ for Archer)'
			update_completed_runs_list(completed_run_list_path,completed_runs)
			exit()

		print "[%s] - launching run_analysis.py ..." % time.strftime("%Y/%m/%d]  [%H:%M")
		completed_runs[folder] = 'run_analysis.py in progress...'
		update_completed_runs_list(completed_run_list_path,completed_runs)
		subprocess.call(['python','%s/run_analysis.py' % pipeline_folder,'--run',run_folder])
		
		print "[%s] - run analysis is completed" % time.strftime("%Y/%m/%d]  [%H:%M")
		completed_runs[folder] = 'completed'
		update_completed_runs_list(completed_run_list_path,completed_runs)

if not new_analysis:
	print "[%s] nothing to do" % time.strftime("%Y/%m/%d]  [%H:%M")
