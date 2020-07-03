#!/usr/bin/python
import re
import os
import subprocess
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-i', '--input-folder',	help="Illumina folder wich contains BCL",	dest='input_folder')
parser.add_option('-o', '--output-folder',	help="Run folder wich will contains FASTQ",	dest='output_folder') 
(options, args) = parser.parse_args()

sbt_pattern = r"^.*[A-Z]{2}[0-9]{3}[A-Z]{1}$" ## ex BESCOND-AW210F (2 lettres 3 chiffres 1 lettre)
sbt_old_pattern = r"^.*[A-Z]{1}[0-9]{3}[A-Z]{1}$" ## ex GUSMINI-Z492F (1 lettre 3 chiffres 1 lettre)
hemato_pattern = r"^.*[0-9]{3}\.[0-9]{3}$" ## ex ETIENNE-MICHELE-204.066 (se termine par 3 chiffres 1 point 3 chiffres)
hemato_old_pattern = r"^.*[0-9]{3}\.[0-9]{2}$" ## ex ETIENNE-MICHELE-204.06 (se termine par 3 chiffres 1 point 2 chiffres)
hemato_retard_pattern = r"^.*[0-9]{3}-[0-9]{3}$" ## ex ETIENNE-MICHELE-204-066 (se termine par 3 chiffres 1 point 3 chiffres) because fuck you
hemato_retard_old_pattern = r"^.*[0-9]{3}-[0-9]{2}$" ## ex ETIENNE-MICHELE-204-06 (se termine par 3 chiffres 1 point 2 chiffres) because fuck you
hemato_long_pattern = r"^.*TU\d\d\.c\d\.R\d\.B[0-9]{3}\.[0-9]{3}$" ## ex LACARRA-FRANCOIS-TU00.c0.R0.B201.069 (se termine par 3 chiffres 1 point 3 chiffres)

# BCL2FASTQ

sample_sheet = '%s/SampleSheet.csv' % option.input_folder
if not os.path.isfile(sample_sheet):
	print "ERROR : No SampleSheet found, aborting"
	exit()

cmd = subprocess.Popen(['/usr/local/bcl2fastq2-v2.2.x/./bin/bcl2fastq',
'--barcode-mismatches','1',
'--mask-short-adapter-reads','0',
'--use-bases-mask','Y150,I8,Y10,Y150',
'-r','2',
'-w','2',
'-p','14',
'-R',options.input_folder,
'--no-lane-splitting'
],
stdout=open('%s/run_bcl2fastq_illumina.stdout.txt' % options.input_folder,'w'),
stderr=open('%s/run_bcl2fastq_illumina.stderr.txt' % options.input_folder,'w')
)
cmd.communicate()

# CREATE OUTPUT FOLDER AND BARCODES.JSON
subprocess.call(['mkdir',option.output_folder]) # pour l'instant

sh_reader = open(sample_sheet,'r')
for line in sh_reader:
	# if line.startswith('Run Name'):
		# run_name = line.split(',')[-1]
		# subprocess.call(['mkdir',option.output_folder])
	while not line.startswith('Sample_ID,'):
		continue
	data = line.split(',')
	sample_id = data[0]
	sample_project = data[3]

	# PARSE SAMPLE NAME ETC
