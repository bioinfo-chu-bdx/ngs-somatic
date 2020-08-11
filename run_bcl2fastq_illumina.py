#!/usr/bin/python
import re
import os
import time
import glob
import json
import uuid
import subprocess
from optparse import OptionParser

# THIS SCRIPT is auto-launched on illumina output folder after sequencing (see Systemd demon)
# USAGE : python run_bcl2fastq_illumina.py /path/to/illumina_output_folder

parser = OptionParser()
parser.add_option('-i', '--illumina-output',help="Illumina folder produced by sequencer wich contains BCL",dest='illumina_folder')
(options, args) = parser.parse_args()

output_location = '/media/n06lbth/sauvegardes_pgm/LAM/panel-capture-v1'

sbt_pattern			= r"^.*[A-Z]{2}[0-9]{3}[A-Z]{1}$" ## ex BESCOND-AW210F (2 lettres 3 chiffres 1 lettre)
sbt_old_pattern		= r"^.*[A-Z]{1}[0-9]{3}[A-Z]{1}$" ## ex GUSMINI-Z492F (1 lettre 3 chiffres 1 lettre)
hemato_pattern		= r"^.*[0-9]{3}\.[0-9]{3}$" ## ex ETIENNE-MICHELE-204.066 (se termine par 3 chiffres 1 point 3 chiffres)
hemato_old_pattern	= r"^.*[0-9]{3}\.[0-9]{2}$" ## ex ETIENNE-MICHELE-204.06 (se termine par 3 chiffres 1 point 2 chiffres)
hemato_pattern2		= r"^.*[0-9]{3}-[0-9]{3}$" ## ex ETIENNE-MICHELE-204-066 (se termine par 3 chiffres 1 tiret 3 chiffres)
hemato_old_pattern2	= r"^.*[0-9]{3}-[0-9]{2}$" ## ex ETIENNE-MICHELE-204-06 (se termine par 3 chiffres 1 tiret 2 chiffres)
hemato_long_pattern	= r"^.*TU\d\d\.c\d\.R\d\.B[0-9]{3}\.[0-9]{3}$" ## ex LACARRA-FRANCOIS-TU00.c0.R0.B201.069 (se termine par 3 chiffres 1 point 3 chiffres)
hemato_long_pattern2= r"^.*[0-9]{3}\.[0-9]{4}$" ## ex ETIENNE-MICHELE-204.0066 (se termine par 3 chiffres 1 point 4 chiffres)
hemato_long_pattern3= r"^.*[0-9]{3}\-[0-9]{4}$" ## ex ETIENNE-MICHELE-204.0066 (se termine par 3 chiffres 1 tiret 4 chiffres)

control_names = ['H2O','H20','NTC','ACROMETRIX','BAF5','BAF-5','HD300','HD301','HD802','HD901','HD748','HORIZON','TEMOIN'] # liste des noms possibles pour les temoins negatifs

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

#### TODO ####
# A FAIRE DANS UN AUTRE SCRIPT, AVEC SYSTEMD :
# trouver comment verifier que nouveau run est apparu
# pour savoir si nouveau, utiliser la BDD SQL voir si deja dedans?
# ensuite, savoir quand nouveau run est termine : est-ce que sampleSheet apparait?


sample_sheet = '%s/SampleSheet.csv' % options.illumina_folder
if not os.path.isfile(sample_sheet):
	print "ERROR : No SampleSheet found, aborting"
	exit()

# BCL2FASTQ
cmd = subprocess.Popen(['/usr/local/bcl2fastq2-v2.2.x/./bin/bcl2fastq',
'--barcode-mismatches','1',
'--mask-short-adapter-reads','0',
'--use-bases-mask','Y150,I8,Y10,Y150',
'-r','2',
'-w','2',
'-p','14',
'-R',options.illumina_folder,
'--no-lane-splitting'
],
stdout=open('%s/run_bcl2fastq_illumina.stdout.txt' % options.illumina_folder,'w'),
stderr=open('%s/run_bcl2fastq_illumina.stderr.txt' % options.illumina_folder,'w'))
cmd.communicate()

# PARSE SAMPLE SHEET
print "- PARSING SAMPLESHEET"
ss_reader = open(sample_sheet,'r')
experiment_name = 'unnamed_run_%s' % time.strftime("%d-%m-%Y-%Hh%M")
sample_section = False
samples = {}
number = 1
for line in ss_reader:
	if line in ['\n','\r\n']:
		continue
	if line.startswith('Experiment Name'):
		experiment_name = line.replace('\r\n','').split(',')[-1]
	if line.startswith('Platform'):
		platform = line.replace('\r\n','').split(',')[-1]
	if line.startswith('Target technique'):
		target_technique = line.replace('\r\n','').split(',')[-1]
	elif line.startswith('Sample_ID,'):
		sample_section = True
		continue
	if sample_section:
		sample_data = line.replace('\r\n','').split(',')
		sample_id = sample_data[0]
		samples[sample_id] = {}
		samples[sample_id]['sample_project'] = sample_data[3]
		samples[sample_id]['sample_number'] = 'S%s' % number
		number+=1

# CREATE OUTPUT FOLDER 
print "- CREATING RUN FOLDER"
run_folder = '%s/%s' % (output_location,experiment_name)
while os.path.isdir(run_folder):
	run_folder = '%s_%s' % (run_folder,time.strftime("%d-%m-%Y-%Hh%M"))
print "\t - %s" % (run_folder)
subprocess.call(['mkdir',run_folder])

# CREATE PATIENT FOLDERS & TRANSFER FASTQ
print "- CREATING PATIENT FOLDERS & TRANSFER FASTQ"
for sample in samples.keys():
	print "\t - %s" % (sample)
	sample_folder = '%s/%s' % (run_folder,sample)
	if not os.path.exists(sample_folder):
		subprocess.call(['mkdir', sample_folder])
	# find all sample fastqs
	sample_fastqs = glob.glob('%s/Data/Intensities/BaseCalls/%s/%s*.fastq.gz' % (options.illumina_folder,samples[sample_id]['sample_project'],sample))
	for fastq in sample_fastqs:
		print "\t\t - %s found" % (os.path.basename(fastq))
		subprocess.call(['cp', fastq, sample_folder])

# AND GENERATE BARCODES.JSON
print "- GENERATE BARCODES.JSON"
barcodes_json = {}
for sample in samples.keys():
	barcode = samples[sample]['sample_number']
	target = global_param['run_type'][samples[sample]['sample_project']]['target_bed'].split('/')[-1]

	id_in_filename = True
	iscontrol = False
	for control_name in control_names:
		if control_name in sample.upper():
			if re.match(sbt_pattern,sample) or re.match(sbt_old_pattern,sample): # control should not have dna number... this is to avoid REVLEG-AH200F
				continue
			else:
				iscontrol = True
	if iscontrol:
		print "\t -is control"
		random_uuid = uuid.uuid1()
		sample_id = 'CONTROL-%s' % random_uuid.hex[:8].upper()
		id_in_filename = False
	if not iscontrol:
		if re.match(sbt_pattern,sample):
			print "\t -sbt pattern ('XX123X')"
			sample_id = sample[-6:]
			sample = sample[:-7]
		elif re.match(sbt_old_pattern,sample):
			print "\t -sbt old pattern ('X123X')"
			sample_id = sample[-5:]
			sample = sample[:-6]
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

	# A VOIR SI BESOIN - ELIMINATE sample END - or _ 
	# A VOIR SI BESOIN - CHECK AND CORRECT SAMPLE ID DUPLICATE
	# A VOIR SI BESOIN - CHECK AND CORRECT SAMPLE DUPLICATE

	if id_in_filename:
		sample_name = '%s_%s' % (sample,sample_id)
	else:
		sample_name = sample
	barcodes_json[barcode] = {'reference':'hg19','sample':sample_name,'sample_id':'%s' % sample_id,'barcode_description':'','target_region_filepath':'%s' % target,'platform':platfom,'target_technique':target_technique}

json_text = json.dumps(barcodes_json, indent=4)
bc_json = open('%s/barcodes.json' % run_folder,'w')
bc_json.write(json_text)
bc_json.close()
