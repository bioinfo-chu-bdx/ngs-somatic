#!/usr/bin/python
import os
import sys
import time
import glob
import json
import zipfile
import logging
import subprocess
from datetime import date
from optparse import OptionParser


def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

def check_stderr(stderr_path):
	with open(stderr_path,'r') as stderr:
		for line in stderr.readlines():
			logging.info(line.replace('\n',''))

def init_log(log_path): # set up log file
	logging.basicConfig(level=logging.INFO,format='%(levelname)-7s %(message)s',filename='/%s/analysis.log' % log_path,filemode='w')
	console = logging.StreamHandler()
	console.setLevel(logging.INFO)
	logging.getLogger('').addHandler(console) # add the handler to the root logger


### GATHERING PARAMETERS ############################################################

FNULL = open(os.devnull, 'w')
parser = OptionParser()
parser.add_option('-r', '--run',				help="run folder path for FULL RUN ANALYSIS",dest='run') 
parser.add_option('-s', '--sample',				help="sample folder path for SINGLE SAMPLE ANALYSIS", dest='sample')
parser.add_option('-x', '--skip-preprocessing',	help="start analysis directly with BAM (skip fastq pre-processing and alignment)",dest='skip_preprocessing',default=False,action='store_true')
parser.add_option('-y', '--skip-calling',		help="start analysis directly with Annotation and Finalreport",dest='skip_calling',default=False,action='store_true')
(options, args) = parser.parse_args()

if options.run and options.sample:
	sys.stderr.write("[run_analysis.py] Error: choose either <--run> or <--sample>, not both\n")
	sys.exit()
if options.run:
	run_folder = options.run
	init_log(run_folder)
elif options.sample:
	if options.sample.endswith('/'):
		options.sample = options.sample[:-1]
	sample_folder = options.sample
	run_folder = os.path.dirname(sample_folder)
	init_log(sample_folder)
else:
	sys.stderr.write("[run_analysis.py] Error: no <--fastq> or <--full-run> specified\n")
	sys.exit()

if os.path.isfile('%s/barcodes.json' % run_folder):
	with open('%s/barcodes.json' % run_folder, 'r') as g:
		barcodes_json = json.load(g)
else:
	logging.error("barcodes.json file not found")
	sys.exit()

if options.sample:
	for barcode in barcodes_json.keys():
		if not barcodes_json[barcode]['sample'] == options.sample.split('/')[-1]:
			del barcodes_json[barcode]

ordered_barcodes = [item[1] for item in sorted([(barcodes_json[barcode]['sample'],barcode) for barcode in barcodes_json])]

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

sampleask = False
######################################################################################

logging.info(" [%s] STARTING ANALYSIS ..." % (time.strftime("%H:%M:%S")))

### RUN ANALYSIS ######################################################################
for barcode in ordered_barcodes:
	barcodes_json[barcode]['fastq'] = '%s/%s/%s_%s_R1_001.fastq.gz' % (run_folder,barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode)
	barcodes_json[barcode]['bam'] = '%s/%s/%s_%s.bam' % (run_folder,barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode)
	barcodes_json[barcode]['sample_folder'] = '%s/%s' % (run_folder,barcodes_json[barcode]['sample'])
	barcodes_json[barcode]['intermediate_folder'] = '%s/intermediate_files' % barcodes_json[barcode]['sample_folder']
	barcodes_json[barcode]['reference'] = global_param['default_reference']
	barcodes_json[barcode]['target_bed'] = global_param['run_type'][barcodes_json[barcode]['project']]['target_bed']
	barcodes_json[barcode]['covered_bed'] = global_param['run_type'][barcodes_json[barcode]['project']].get('covered_bed','target_bed') # si covered existe, sinon target
	barcodes_json[barcode]['merged_bed'] = global_param['run_type'][barcodes_json[barcode]['project']]['merged_bed']
	barcodes_json[barcode]['intervals'] = global_param['run_type'][barcodes_json[barcode]['project']].get('intervals',False)
	if barcodes_json[barcode]['covered_bed'] == barcodes_json[barcode]['target_bed']:
		logging.info("WARNING : covered_bed is also roi_bed for %s" % barcodes_json[barcode]['sample'])
	if not barcodes_json[barcode]['intervals']:
		logging.info("WARNING : intervals file missing for %s" % barcodes_json[barcode]['sample'])

	# DOSSIER INTERMEDIAIRE
	if os.path.isfile('%s.zip' % barcodes_json[barcode]['intermediate_folder']): # UNZIP ALL intermediates_files.zip (if RERUN)
		fnull = open(os.devnull, 'w')
		subprocess.call(['unzip', '-o', '%s.zip' % barcodes_json[barcode]['intermediate_folder'],'-d', barcodes_json[barcode]['intermediate_folder']],stdout=fnull)
	elif not os.path.isdir(barcodes_json[barcode]['intermediate_folder']):
		subprocess.call(['mkdir', barcodes_json[barcode]['intermediate_folder']])

	logging.info("\t - %s, analysisID = %s, panel = %s" % (barcodes_json[barcode]['sample'],barcodes_json[barcode]['analysis_id'],barcodes_json[barcode]['project']))

###   __   __   ___     __   __   __   __   ___  __   __          __    ###
###  |__) |__) |__  __ |__) |__) /  \ /  ` |__  /__` /__` | |\ | / _`   ###
###  |    |  \ |___    |    |  \ \__/ \__, |___ .__/ .__/ | | \| \__>   ###

logging.info("\n- [%s] PRE-PROCESSING ..." % (time.strftime("%H:%M:%S")))

if not options.skip_preprocessing:
	for barcode in ordered_barcodes:
		logging.info("\t- %s :" % barcodes_json[barcode]['sample'])
		prep_folder = '%s/pre-processing' % barcodes_json[barcode]['intermediate_folder']
		if not os.path.isdir(prep_folder):
			subprocess.call(['mkdir', prep_folder])

		paired_end_classic = False
		paired_end_mbc = False
		# RUN PAIRED END WITH MBC
		if os.path.isfile(barcodes_json[barcode]['fastq'].replace('R1_001','R3_001')):
			paired_end_mbc = True
			barcodes_json[barcode]['fastq_r2'] = barcodes_json[barcode]['fastq'].replace('R1_001','R2_001')
			barcodes_json[barcode]['fastq_r3'] = barcodes_json[barcode]['fastq'].replace('R1_001','R3_001')
		# RUN PAIRED END NO MBC
		elif os.path.isfile(barcodes_json[barcode]['fastq'].replace('R1_001','R2_001')):
			paired_end_classic = True
			barcodes_json[barcode]['fastq_r2'] = barcodes_json[barcode]['fastq'].replace('R1_001','R2_001')
		barcodes_json[barcode]['fastq_r1'] = barcodes_json[barcode]['fastq']

		# (OPTIONAL) AGENT TRIM
		if paired_end_mbc:
			logging.info("\t\t- [%s] Trim MBC ..." % time.strftime("%H:%M:%S"))
			cmd = subprocess.Popen(['%s/agent/agent.sh' % pipeline_folder,'trim','-xt','-fq1',barcodes_json[barcode]['fastq_r1'],'-fq2',barcodes_json[barcode]['fastq_r3'],'-out_loc',prep_folder], stdout=open('%s/trim.stdout.txt' % prep_folder,'w'), stderr=open('%s/trim.stderr.txt' % prep_folder,'w'))
			cmd.communicate()
			barcodes_json[barcode]['fastq_r1'] = glob.glob('%s/*_R1_*fastq*' % prep_folder)[0]
			barcodes_json[barcode]['fastq_r3'] = glob.glob('%s/*_R3_*fastq*' % prep_folder)[0]

		# ALIGNMENT - FASTQ TO SAM
		logging.info("\t\t- [%s] BWA-MEM Alignment ..." % time.strftime("%H:%M:%S"))
		read_group = '@RG\\tID:%s\\tLB:%s\\tPL:%s\\tPU:%s\\tSM:%s' % ('%s' % barcodes_json[barcode]['sample_id'],'Target-Myeloid-v1','illumina','NDX550372',barcodes_json[barcode]['sample']) # READ GROUP (which reflects which library a read belongs to and what lane it was sequenced in on the flowcell)
		if paired_end_mbc :
			cmd = subprocess.Popen(['%s/bwa/bwa' % pipeline_folder,'mem','-C','-t','12','-M','-R',read_group,'-v','1','-o','%s/%s.sam' % (prep_folder,barcodes_json[barcode]['sample']),barcodes_json[barcode]['reference'],barcodes_json[barcode]['fastq_r1'],barcodes_json[barcode]['fastq_r3']], stdout=open('%s/bwa_mem.stdout.txt' % prep_folder,'w'), stderr=open('%s/bwa_mem.stderr.txt' % prep_folder,'w'))
		elif paired_end_classic:
			cmd = subprocess.Popen(['%s/bwa/bwa' % pipeline_folder,'mem','-C','-t','12','-M','-R',read_group,'-v','1','-o','%s/%s.sam' % (prep_folder,barcodes_json[barcode]['sample']),barcodes_json[barcode]['reference'],barcodes_json[barcode]['fastq_r1'],barcodes_json[barcode]['fastq_r2']], stdout=open('%s/bwa_mem.stdout.txt' % prep_folder,'w'), stderr=open('%s/bwa_mem.stderr.txt' % prep_folder,'w'))
		else:
			cmd = subprocess.Popen(['%s/bwa/bwa' % pipeline_folder,'mem','-C','-t','12','-M','-R',read_group,'-v','1','-o','%s/%s.sam' % (prep_folder,barcodes_json[barcode]['sample']),barcodes_json[barcode]['reference'],barcodes_json[barcode]['fastq_r1']], stdout=open('%s/bwa_mem.stdout.txt' % prep_folder,'w'), stderr=open('%s/bwa_mem.stderr.txt' % prep_folder,'w'))
		cmd.communicate()

		if paired_end_classic :
			## MARKDUPLICATESSPARK (& SORT)
			logging.info("\t\t- [%s] Marking duplicates and sorting BAM ..." % time.strftime("%H:%M:%S"))
			cmd = subprocess.Popen(['docker','exec','-it','gatk','gatk','MarkDuplicatesSpark','--verbosity','ERROR','-I','%s/%s.sam' % (prep_folder,barcodes_json[barcode]['sample']),'-O','%s/%s.markdup.sorted.bam' % (prep_folder,barcodes_json[barcode]['sample'])], stdout=open('%s/markduplicates.stdout.txt' % prep_folder,'w'), stderr=open('%s/markduplicates.stderr.txt' % prep_folder,'w'))
			cmd.communicate()

			# QUALITY SCORE RECALIBRATION (KEEP KNOWN SITES?)
			logging.info("\t\t- [%s] Base recalibration ..." % time.strftime("%H:%M:%S"))
			cmd = subprocess.Popen(['docker','exec','-it','gatk','gatk','BaseRecalibratorSpark','-I','%s/%s.markdup.sorted.bam' % (prep_folder,barcodes_json[barcode]['sample']),'-R',barcodes_json[barcode]['reference'],'--known-sites','%s/reference_files/mutect/af-only-gnomad.raw.sites.b37.vcf.gz' % pipeline_folder,'-O','%s/%s.recal_data.table' % (prep_folder,barcodes_json[barcode]['sample'])], stdout=open('%s/baserecalibrator.stdout.txt' % prep_folder,'w'), stderr=open('%s/baserecalibrator.stderr.txt' % prep_folder,'w'))
			cmd.communicate()
			cmd = subprocess.Popen(['docker','exec','-it','gatk','gatk','ApplyBQSRSpark','-I','%s/%s.markdup.sorted.bam' % (prep_folder,barcodes_json[barcode]['sample']),'-R',barcodes_json[barcode]['reference'],'--bqsr-recal-file','%s/%s.recal_data.table' % (prep_folder,barcodes_json[barcode]['sample']),'-O',barcodes_json[barcode]['bam']], stdout=open('%s/applybqsr.stdout.txt' % prep_folder,'w'), stderr=open('%s/applybqsr.stderr.txt' % prep_folder,'w'))
			cmd.communicate()
		else:
			if paired_end_mbc:
				# PROCESS MBC
				logging.info("\t\t- [%s] Process MBC ..." % time.strftime("%H:%M:%S"))
				cmd = subprocess.Popen(['%s/agent/agent.sh' % pipeline_folder,'locatit','-i','-C','-m','1','-U','-r','-IS','-OB','-c','100','-q','20','-c','100','-l',barcodes_json[barcode]['covered_bed'],'-o','%s/%s.markdup.sorted.bam' % (prep_folder,barcodes_json[barcode]['sample']),'%s/%s.sam' % (prep_folder,barcodes_json[barcode]['sample']),barcodes_json[barcode]['fastq_r2']], stdout=open('%s/locatit.stdout.txt' % prep_folder,'w'), stderr=open('%s/locatit.stderr.txt' % prep_folder,'w'))
				cmd.communicate()
			else:
				# SAM TO BAM
				logging.info("\t\t- [%s] Converting SAM to BAM ..." % time.strftime("%H:%M:%S"))
				cmd = subprocess.Popen(['samtools','view','-@','12','-bS','-o','%s/%s.markdup.sorted.bam' % (prep_folder,barcodes_json[barcode]['sample']),'%s/%s.sam' % (prep_folder,barcodes_json[barcode]['sample'])], stdout=open('%s/samtools_view.stdout.txt' % prep_folder,'w'), stderr=open('%s/samtools_view.stderr.txt' % prep_folder,'w'))
				cmd.communicate()

			# SORT BAM
			logging.info("\t\t- [%s] Sorting BAM ..." % time.strftime("%H:%M:%S"))
			cmd = subprocess.Popen(['samtools','sort','-@','12','-O','bam','-o',barcodes_json[barcode]['bam'],'%s/%s.markdup.sorted.bam' % (prep_folder,barcodes_json[barcode]['sample'])], stdout=open('%s/samtools_sort.stdout.txt' % prep_folder,'w'), stderr=open('%s/samtools_sort.stderr.txt' % prep_folder,'w'))
			cmd.communicate()

			# INDEX
			logging.info("\t\t- [%s] Indexing BAM ..." % time.strftime("%H:%M:%S"))
			cmd = subprocess.Popen(['samtools','index',barcodes_json[barcode]['bam']], stdout=open('%s/samtools_index.stdout.txt' % prep_folder,'w'), stderr=open('%s/samtools_index.stderr.txt' % prep_folder,'w'))
			cmd.communicate()

		# remove SAM and BAM temp files (if bam correctly produced)
		if os.path.exists(barcodes_json[barcode]['bam']):
			subprocess.call(['rm','%s/%s.sam' % (prep_folder,barcodes_json[barcode]['sample'])],stdout=FNULL,stderr=FNULL)
			subprocess.call(['rm','%s/%s.markdup.sorted.bam' % (prep_folder,barcodes_json[barcode]['sample'])],stdout=FNULL,stderr=FNULL)
			subprocess.call(['rm','%s/%s.markdup.sorted.bam.bai' % (prep_folder,barcodes_json[barcode]['sample'])],stdout=FNULL,stderr=FNULL)
			subprocess.call(['rm','%s/%s.markdup.sorted.bam.sbi' % (prep_folder,barcodes_json[barcode]['sample'])],stdout=FNULL,stderr=FNULL)


##### #####            __              ___     __                         __                __                     __  ___      ___    __       ##### ##### 
##### ##### \  /  /\  |__) |  /\  |\ |  |  __ /  `  /\  |    |    | |\ | / _`     /\  |\ | |  \     /\  |\ | |\ | /  \  |   /\   |  | /  \ |\ | ##### ##### 
##### #####  \/  /~~\ |  \ | /~~\ | \|  |     \__, /~~\ |___ |___ | | \| \__>    /~~\ | \| |__/    /~~\ | \| | \| \__/  |  /~~\  |  | \__/ | \| ##### ##### 
##### #####                                                                                                                                     ##### ##### 

logging.info("\n- [%s] VARIANT-CALLING AND ANNOTATION ..." % (time.strftime("%H:%M:%S")))
for barcode in ordered_barcodes:
	logging.info("\t- %s :" % (barcodes_json[barcode]['sample']))

		###            __   __     __  ___ 
		### \  /  /\  |__) |  \ | /  `  |  
		###  \/  /~~\ |  \ |__/ | \__,  |  

	logging.info("\t\t - [%s] VarDict ..." % time.strftime("%H:%M:%S"))
	vardict_folder = '%s/vardict' % barcodes_json[barcode]['intermediate_folder']
	if not os.path.isdir(vardict_folder):
		subprocess.call(['mkdir', vardict_folder])

	vardict_chunk = 12
	logging.info("\t\t\t - [%s] split_bed.py ..." % time.strftime("%H:%M:%S"))
	os.chdir(vardict_folder)
	cmd = subprocess.Popen(['python','%s/scripts/split_bed.py' % pipeline_folder,'--bed',barcodes_json[barcode]['target_bed'],'--scatter-count',str(vardict_chunk),'--output-folder',vardict_folder],stdout=open('%s/split_bed.stdout.txt'%vardict_folder,'w'),stderr=open('%s/split_bed.stderr.txt'%vardict_folder,'w'))
	cmd.communicate()

	ps = []
	vcf_chunk_list = []
	logging.info("\t\t\t - [%s] Running %s chunks ..." % (time.strftime("%H:%M:%S"),vardict_chunk))
	FNULL = open(os.devnull, 'w')
	for i in range(vardict_chunk):
		bed_chunk = '%s/%04d-scattered.bed'%(vardict_folder,i)
		vcf_chunk = '%s/%s.vardict.%04d.vcf' % (vardict_folder,barcodes_json[barcode]['sample'],i)
		vcf_chunk_list.append(vcf_chunk)
		cmd1 = subprocess.Popen(['perl','%s/vardict/VarDictJava/VarDict/vardict' % pipeline_folder,
			'-G', barcodes_json[barcode]['reference'],
			'-f','0.01',
			'-N',barcodes_json[barcode]['sample'],
			'-b', barcodes_json[barcode]['bam'],
			'-c','1','-S','2','-E','3','-g','4',
			bed_chunk],
			stdout=open('%s/vardict.%04d.stdout.txt' % (vardict_folder,i),'w'),
			stderr=open('%s/vardict.%04d.stderr.txt' % (vardict_folder,i),'w'))
		ps.append(cmd1)
	for cmd1 in ps:
		cmd1.wait()

	logging.info("\t\t\t - [%s] teststrandbias.R ..." % time.strftime("%H:%M:%S"))
	ps = []
	for i in range(vardict_chunk):
		cmd2 = subprocess.Popen(['%s/vardict/VarDictJava/VarDict/teststrandbias.R' % pipeline_folder],
		stdin=open('%s/vardict.%04d.stdout.txt' % (vardict_folder,i),'r'),
		stdout=open('%s/teststrandbias.%04d.stdout.txt' % (vardict_folder,i),'w'),
		stderr=open('%s/teststrandbias.%04d.stderr.txt' % (vardict_folder,i),'w'))
		ps.append(cmd2)
	for cmd2 in ps:
		cmd2.wait()

	logging.info("\t\t\t - [%s] var2vcf_valid.pl ..." % time.strftime("%H:%M:%S"))
	ps = []
	for i in range(vardict_chunk):
		cmd3 = subprocess.Popen(['perl','%s/vardict/VarDictJava/VarDict/var2vcf_valid.pl' % pipeline_folder,
			'-N',barcodes_json[barcode]['sample'],
			'-E',
			'-f','0.01'],
		stdin=open('%s/teststrandbias.%04d.stdout.txt' % (vardict_folder,i),'r'),
		stdout=open('%s/%s.vardict.%04d.vcf' % (vardict_folder,barcodes_json[barcode]['sample'],i),'w'),
		stderr=open('%s/var2vcf_valid.%04d.stderr.txt' % (vardict_folder,i),'w'))
		ps.append(cmd3)
	for cmd3 in ps:
		cmd3.wait()

	logging.info("\t\t\t - [%s] UpdateVCFSequenceDictionary ..." % time.strftime("%H:%M:%S"))
	ps = []
	for vcf_chunk in vcf_chunk_list:
		cmd = subprocess.Popen(['docker','exec','-it','gatk','gatk','UpdateVCFSequenceDictionary',
			'-V',vcf_chunk,
			'--source-dictionary',barcodes_json[barcode]['bam'],
			'--output',vcf_chunk.replace('.vcf','.contig.vcf')],
		stdout=open('%s/UpdateVCFSequenceDictionary.stdout.txt' % vardict_folder,'w'),
		stderr=open('%s/UpdateVCFSequenceDictionary.stderr.txt' % vardict_folder,'w'))
		ps.append(cmd)
	for cmd in ps:
		cmd.wait()

	logging.info("\t\t\t - [%s] GatherVcfs ..." % time.strftime("%H:%M:%S"))
	args = ['docker','exec','-it','gatk','gatk','GatherVcfs','-O','%s/%s.vardict.vcf' % (vardict_folder,barcodes_json[barcode]['sample'])]
	for vcf_chunk in vcf_chunk_list:
		args.append('-I')
		args.append(vcf_chunk.replace('.vcf','.contig.vcf'))
	cmd = subprocess.Popen(args,stdout=open('%s/gathervcf.stdout.txt' % vardict_folder,'w'),stderr=open('%s/gathervcf.stderr.txt' % vardict_folder,'w'))
	cmd.wait()
	subprocess.call(['stty','sane'])
