#!/usr/bin/python
import re
import os
import time
import glob
import json
import uuid
import subprocess
from optparse import OptionParser

# THIS SCRIPT should be launched on illumina output folder after sequencing (see Systemd demon for automation)
# FROM BCL -> RUN FOLDER WITH FASTQ
# USAGE : python run_bcl2fastq_illumina.py /path/to/illumina_output_folder

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

parser = OptionParser()
parser.add_option('-i', '--illumina_folder',help="Path to illumina folder containing BCL. It is the folder produced by the sequencer",dest='illumina_folder')
(options, args) = parser.parse_args()

sbt_pattern			= r"^.*[A-Z]{2}[0-9]{3}[A-Z]{1}$" ## ex BESCOND-AW210F (2 lettres 3 chiffres 1 lettre)
sbt_old_pattern		= r"^.*[A-Z]{1}[0-9]{3}[A-Z]{1}$" ## ex GUSMINI-Z492F (1 lettre 3 chiffres 1 lettre)
sbt_cdna_pattern1 = r"^.*[A-Z]{3}[0-9]{4}$" ## ex ALIE-PLR1146 (3 lettres 4 chiffres)
sbt_cdna_pattern2 = r"^.*[A-Z]{3}[0-9]{3}$" ## ex ALLARD-PLR905 (3 lettres 3 chiffres)
sbt_cdna_pattern3 = r"^.*[A-Z]{2}[0-9]{2}[A-Z]{1}[0-9]{1}$" ## ex VANDENBERGHE-PL05C2 (2 lettres 2 chiffres 1 lettre 1 chiffre)
hemato_pattern		= r"^.*[0-9]{3}\.[0-9]{3}$" ## ex ETIENNE-MICHELE-204.066 (se termine par 3 chiffres 1 point 3 chiffres)
hemato_old_pattern	= r"^.*[0-9]{3}\.[0-9]{2}$" ## ex ETIENNE-MICHELE-204.06 (se termine par 3 chiffres 1 point 2 chiffres)
hemato_pattern2		= r"^.*[0-9]{3}-[0-9]{3}$" ## ex ETIENNE-MICHELE-204-066 (se termine par 3 chiffres 1 tiret 3 chiffres)
hemato_old_pattern2	= r"^.*[0-9]{3}-[0-9]{2}$" ## ex ETIENNE-MICHELE-204-06 (se termine par 3 chiffres 1 tiret 2 chiffres)
hemato_long_pattern	= r"^.*TU\d\d\.c\d\.R\d\.B[0-9]{3}\.[0-9]{3}$" ## ex LACARRA-FRANCOIS-TU00.c0.R0.B201.069 (se termine par 3 chiffres 1 point 3 chiffres)
hemato_long_pattern2= r"^.*[0-9]{3}\.[0-9]{4}$" ## ex ETIENNE-MICHELE-204.0066 (se termine par 3 chiffres 1 point 4 chiffres)
hemato_long_pattern3= r"^.*[0-9]{3}\-[0-9]{4}$" ## ex ETIENNE-MICHELE-204.0066 (se termine par 3 chiffres 1 tiret 4 chiffres)
hemato_glims_pattern= r"^.*HADN[0-9]{3}\-[0-9]{4}$" ## ex ETIENNE-MICHELE-HADN204-0066
hemato_glims_pattern2= r"^.*HADN[0-9]{3}\.[0-9]{4}$" ## ex ETIENNE-MICHELE-HADN204.0066

# liste des noms possibles pour les temoins negatifs
control_names = ['H2O','H20','NTC','ACROMETRIX','BAF5','BAF-5','HD300','HD301','HD802','HD901','HD748','HORIZON','TEMOIN']

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

ngs_folder = global_param['ngs_results_folder']

sample_sheet = '%s/SampleSheet.csv' % options.illumina_folder
if not os.path.isfile(sample_sheet):
	print "ERROR : No SampleSheet found, aborting"
	exit()

# # PARSE SAMPLE SHEET
print "- PARSING SAMPLESHEET"
run_project = ''
sub_project = ''
ss_reader = open(sample_sheet,'r')
experiment_name = 'unnamed_run_%s' % time.strftime("%d-%m-%Y-%Hh%M")
sample_section = False
samples = {}
number = 1
for line in ss_reader:
	if line in ['\n','\r\n']:
		continue
	if line.startswith('Experiment Name'):
		experiment_name = line.replace('\n','').replace('\r','').split(',')[-1]
	if line.startswith('Platform'):
		platform = line.replace('\n','').replace('\r','').split(',')[-1]
	if line.startswith('System'):
		system = line.replace('\n','').replace('\r','').split(',')[-1]
	if line.startswith('Target technique'):
		target_technique = line.replace('\n','').replace('\r','').split(',')[-1]
	if line.startswith('Run Project'):
		run_project = line.replace('\n','').replace('\r','').split(',')[-1]
	if line.startswith('Sub Project'):
		sub_project = line.replace('\n','').replace('\r','').split(',')[-1]
	elif line.startswith('Sample_ID,'):
		sample_section = True
		continue
	if sample_section:
		sample_data = line.replace('\n','').replace('\r','').split(',')
		sample_id = sample_data[0]
		samples[sample_id] = {}
		samples[sample_id]['panel_name'] = sample_data[3]
		samples[sample_id]['sample_number'] = 'S%s' % number
		number+=1

# FILL TRACABILITE TECHNIQUE DES RUNS CAPTURE
if run_project != 'FusionPlex_CTL':
	try:
		ttr = '/media/n06lbth/sauvegardes_pgm/LAM/EN_LAB_20_2579_Tracabilite_technique_des_runs_capture.xlsx'
		subprocess.call(['python','%s/scripts/tracabilite_technique_capture.py' % pipeline_folder,options.illumina_folder,ttr,experiment_name])
	except:
		pass

# BCL2FASTQ
print "- BCL2FASTQ"
if run_project == 'FusionPlex_CTL':
	bases_mask = 'Y151,I8,I8,Y151'
	print "\t - run FusionPlex_CTL, using bases_mask Y151,I8,I8,Y151"
else:
	bases_mask = 'Y150,I8,Y10,Y150'
	print "\t - run Capture, using bases_mask Y150,I8,Y10,Y150"
cmd = subprocess.Popen(['/usr/local/bcl2fastq2-v2.2.x/./bin/bcl2fastq',
'--barcode-mismatches','1',
'--mask-short-adapter-reads','0',
'-R',options.illumina_folder,
'-r','2', # --loading-threads  : Number of threads to load BCL data. 
'-w','2', # --writing-threads : Number of threads to write FASTQ data. This number must be lower than number of samples
'-p','24', # --processing-threads : Number of threads to process demultiplexing data.
'--use-bases-mask',bases_mask,
'--no-lane-splitting'
],
stdout=open('%s/run_bcl2fastq_illumina.stdout.txt' % options.illumina_folder,'w'),
stderr=open('%s/run_bcl2fastq_illumina.stderr.txt' % options.illumina_folder,'w'))
cmd.communicate()
# '--use-bases-mask','Y150,I8,Y10,Y150',
# When threading is supported, the software uses the follow defaults to manage threads for processing:
#  Four threads for reading the data.
#  Four threads for writing the data.
#  Twenty percent of threads for demultiplexing data.
#  One hundred percent of threads for processing demultiplexed data

# POUR TSO 500 si besoin semble etre Y101;I8;I8;Y101 

# CREATE OUTPUT FOLDER 
print "- CREATING FOLDERS"
output_location = '%s/%s' % (ngs_folder,run_project) # PROJECT FOLDER 
if not os.path.exists(output_location):
	subprocess.call(['mkdir',output_location])
if not sub_project == '': # SUB-PROJECT FOLDER 
	output_location = '%s/%s/%s' % (ngs_folder,run_project,sub_project)
	if not os.path.exists(output_location):
		subprocess.call(['mkdir',output_location])

run_folder = '%s/%s' % (output_location,experiment_name)
# while os.path.isdir(run_folder): ## ?? est-ce vraiment utile
	# run_folder = '%s_%s' % (run_folder,time.strftime("%d-%m-%Y-%Hh%M"))
print "\t - %s" % (run_folder)
subprocess.call(['mkdir',run_folder]) # RUN FOLDER

# CREATE PATIENT FOLDERS & TRANSFER FASTQ
print "- CREATING PATIENT FOLDERS & TRANSFER FASTQ"
if run_project == 'FusionPlex_CTL':
	fastqs = glob.glob('%s/Data/Intensities/BaseCalls/*.fastq.gz' % options.illumina_folder)
	for fastq in fastqs:
		subprocess.call(['mv', fastq, run_folder])
else:
	for sample in samples.keys():
		print "\t - %s" % (sample)
		sample_folder = '%s/%s' % (run_folder,sample)
		if not os.path.exists(sample_folder):
			subprocess.call(['mkdir', sample_folder])
		# find all sample fastqs
		# sample_fastqs = glob.glob('%s/Data/Intensities/BaseCalls/%s/%s*.fastq.gz' % (options.illumina_folder,samples[sample_id]['panel'],sample)) # SI project dans DATA
		sample_fastqs = glob.glob('%s/Data/Intensities/BaseCalls/%s*.fastq.gz' % (options.illumina_folder,sample))
		for fastq in sample_fastqs:
			print "\t\t - %s found" % (os.path.basename(fastq))
			if fastq.startswith('Undetermined_S0'):
				# sample_fastqs.remove(fastq)
				# continue
				subprocess.call(['mv', fastq, run_folder])
			else:
				subprocess.call(['mv', fastq, sample_folder])

# AND GENERATE BARCODES.JSON
print "- GENERATE BARCODES.JSON"
barcodes_json = {}
for sample in samples.keys():
	barcode = samples[sample]['sample_number']
	panel = samples[sample]['panel_name']
	target = global_param['panel'][panel]['target_bed'].split('/')[-1]
	id_in_filename = True
	is_control = 0
	for control_name in control_names:
		if control_name in sample.upper():
			# control should not have dna number... this test is to avoid case like : REVLEG-AH200F
			if re.match(sbt_pattern,sample) or re.match(sbt_old_pattern,sample):
				continue
			else:
				is_control = 1
	if is_control == 1:
		print "\t -is control"
		random_uuid = uuid.uuid1()
		sample_id = 'CONTROL-%s' % random_uuid.hex[:8].upper()
		id_in_filename = False
	else:
		if re.match(hemato_glims_pattern,sample):
			print "\t -hemato glims pattern ('HADN123-4567')"
			sample_id = sample[-12:]
			sample = sample[:-13]
		elif re.match(hemato_glims_pattern2,sample):
			print "\t -hemato glims pattern ('HADN123.4567')"
			sample_id = sample[-12:]
			sample = sample[:-13]
		elif re.match(sbt_pattern,sample):
			print "\t -sbt pattern ('XX123X')"
			sample_id = sample[-6:]
			sample = sample[:-7]
		elif re.match(sbt_old_pattern,sample):
			print "\t -sbt old pattern ('X123X')"
			sample_id = sample[-5:]
			sample = sample[:-6]
		elif re.match(sbt_cdna_pattern1,sample):
			print "\t -sbt cdna pattern1 ('XXX1234')"
			sample_id = sample[-7:]
			sample = sample[:-8]
		elif re.match(sbt_cdna_pattern2,sample):
			print "\t -sbt cdna pattern2 ('XXX123')"
			sample_id = sample[-6:]
			sample = sample[:-7]
		elif re.match(sbt_cdna_pattern3,sample):
			print "\t -sbt cdna pattern3 ('XX12X3')"
			sample_id = sample[-6:]
			sample = sample[:-7]
		elif re.match(hemato_long_pattern,sample):
			print "\t -hemato pattern ('.c0.R0.B123.456')"
			sample_id = sample[-7:]
			sample = sample[:-20]
		elif re.match(hemato_pattern,sample):
			print "\t -hemato pattern ('123.456')"
			sample_id = sample[-7:]
			sample = sample[:-8]
		elif re.match(hemato_old_pattern,sample):
			print "\t -hemato pattern ('123.45')"
			sample_id = sample[-6:]
			sample = sample[:-7]
		elif re.match(hemato_pattern2,sample):
			print "\t -hemato pattern ('123-456')"
			sample_id = sample[-7:]#.replace('-','.')
			sample = sample[:-8]
		elif re.match(hemato_old_pattern2,sample):
			print "\t -hemato pattern ('123-45')"
			sample_id = sample[-6:]#.replace('-','.')
			sample = sample[:-7]
		elif re.match(hemato_long_pattern2,sample):
			print "\t -hemato pattern ('123.4567')"
			sample_id = sample[-8:]
			sample = sample[:-9]
		elif re.match(hemato_long_pattern3,sample):
			print "\t -hemato pattern ('123-4567')"
			sample_id = sample[-8:]#.replace('-','.')
			sample = sample[:-9]
		else: # NO REGEX MATCH
			print "\t -no regex match"
			random_uuid = uuid.uuid1()
			random_uuid = 'ID'+random_uuid.hex[:8].upper()
			sample_id = random_uuid
			id_in_filename = False

	# IF DNA NUMBER = 000.000 -< unknow; generate uuid instead
	if sample_id == '000.000':
		random_uuid = uuid.uuid1()
		random_uuid = 'ID'+random_uuid.hex[:8].upper()
		sample_id = random_uuid
		id_in_filename = False

	while sample.endswith('-') or sample.endswith('_'):
		sample = sample[:-1]

	if id_in_filename:
		sample_full_name = '%s_%s' % (sample,sample_id)
	else:
		sample_full_name = sample

	barcodes_json[barcode] = {'reference':'hg19','sample':sample_full_name,'sample_id':'%s' % sample_id,'is_control':is_control,'library':'','description':'','target_bed':'%s' % target,'platform':platform,'system':system,'target_technique':target_technique,'panel':panel,'analysis_id':''}

json_text = json.dumps(barcodes_json, indent=4, sort_keys=True)
bc_json = open('%s/barcodes.json' % run_folder,'w')
bc_json.write(json_text)
bc_json.close()
