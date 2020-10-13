#!/usr/bin/python
import os
import sys
import time
import glob
import json
import uuid
import sqlite3
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

def normalize_folder_path(folder_path):
	if os.path.isdir(folder_path):
		if folder_path.endswith('/'):
			folder_path = folder_path[:-1]
	return folder_path

def check_stderr(stderr_path):
	with open(stderr_path,'r') as stderr:
		for line in stderr.readlines():
			logging.info(line.replace('\n',''))

def init_log(log_path): # set up log file
	logging.basicConfig(level=logging.INFO,format='%(levelname)-7s %(message)s',filename='%s/analysis.log' % log_path,filemode='a') # /%s/analysis.log?
	console = logging.StreamHandler()
	console.setLevel(logging.INFO)
	logging.getLogger('').addHandler(console) # add the handler to the root logger


### GATHERING PARAMETERS ############################################################

FNULL = open(os.devnull, 'w')
parser = OptionParser()
parser.add_option('-r', '--run',					help="run folder path for FULL RUN ANALYSIS",dest='run') 
parser.add_option('-s', '--sample',					help="sample folder path for SINGLE SAMPLE ANALYSIS", dest='sample')
parser.add_option('-p', '--skip-pre-processing',	help="starts analysis directly with variant calling on BAM (skip fastq pre-processing and alignment)",dest='skip_preprocessing',default=False,action='store_true')
parser.add_option('-c', '--skip-coverage-analysis',	help="skip coverage analysis step",dest='skip_coverage_analysis',default=False,action='store_true')
parser.add_option('-v', '--skip-calling',			help="starts analysis directly with Annotation (skip variant calling and and previous steps)",dest='skip_calling',default=False,action='store_true')
parser.add_option('-i', '--skip-insert-db',			help="s.....s steps)",dest='skip_insert_db',default=False,action='store_true')
parser.add_option('-a', '--skip-annotation',		help="starts analysis directly with Finalreport (skip annotation and previous steps)",dest='skip_annotation',default=False,action='store_true')
(options, args) = parser.parse_args()

if options.skip_annotation:
	options.skip_preprocessing = True
	options.skip_calling = True
elif options.skip_calling:
	options.skip_preprocessing = True

if options.run and options.sample:
	sys.stderr.write("[run_analysis.py] Error: choose either <--run> or <--sample>, not both\n")
	sys.exit()
if options.run:
	options.run = normalize_folder_path(options.run)
	run_folder = options.run
	init_log(run_folder)
elif options.sample:
	options.sample = normalize_folder_path(options.sample)
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

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

######################################################################################

logging.info(" [%s] STARTING ANALYSIS ..." % (time.strftime("%H:%M:%S")))

### RUN ANALYSIS ######################################################################
run_name = os.path.basename(run_folder)
run_date = date.today()
print "- Create new Run DB entry (if inexistant)..."
platform = barcodes_json[list(barcodes_json)[0]]['platform'].lower()
db_cur.execute("INSERT OR IGNORE INTO Run (runID, platform, runPath, runDate) VALUES ('%s', '%s', '%s', '%s')" % (run_name,platform,run_folder,run_date))
db_con.commit()

for barcode in ordered_barcodes:
	barcodes_json[barcode]['fastq'] = '%s/%s/%s_%s_R1_001.fastq.gz' % (run_folder,barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode)
	barcodes_json[barcode]['bam'] = '%s/%s/%s_%s.bam' % (run_folder,barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode)
	barcodes_json[barcode]['sample_folder'] = '%s/%s' % (run_folder,barcodes_json[barcode]['sample'])
	barcodes_json[barcode]['intermediate_folder'] = '%s/intermediate_files' % barcodes_json[barcode]['sample_folder']
	barcodes_json[barcode]['reference'] = global_param['default_reference']
	barcodes_json[barcode]['target_bed'] = global_param['panel'][barcodes_json[barcode]['panel']]['target_bed']
	barcodes_json[barcode]['covered_bed'] = global_param['panel'][barcodes_json[barcode]['panel']].get('covered_bed','target_bed') # si covered existe, sinon target
	barcodes_json[barcode]['merged_bed'] = global_param['panel'][barcodes_json[barcode]['panel']]['merged_bed']
	barcodes_json[barcode]['intervals'] = global_param['panel'][barcodes_json[barcode]['panel']].get('intervals',False)
	barcodes_json[barcode]['tvc_param'] = global_param['panel'][barcodes_json[barcode]['panel']].get('vc_parameters',False)
	barcodes_json[barcode]['tvc_param_hotspot'] = global_param['panel'][barcodes_json[barcode]['panel']].get('vc_parameters_hotspot_only',False)
	barcodes_json[barcode]['hotspot_vcf'] = global_param['panel'][barcodes_json[barcode]['panel']].get('hotspot_vcf',False)
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

	## CREATE SAMPLE DB ENTRY - IF INEXISTANT -
	db_cur.execute("INSERT OR IGNORE INTO Sample (sampleID, sampleName, isControl) VALUES ('%s', '%s', %s)" % (barcodes_json[barcode]['sample_id'],barcodes_json[barcode]['sample'],barcodes_json[barcode]['is_control']))
	db_con.commit()

	## CREATE ANALYSIS DB ENTRY - IF INEXISTANT -
	json_modified = False
	if barcodes_json[barcode]['analysis_id'] == '': # not defined yet
	# verify if not already in DB :
		db_cur.execute("SELECT analysisID FROM Analysis WHERE sample='%s' AND barcode='%s' AND run='%s'" % (barcodes_json[barcode]['sample_id'],barcode,run_name))
		db_analysis = db_cur.fetchone()
		if db_analysis is not None:
			barcodes_json[barcode]['analysis_id'] = db_analysis['analysisID']
			continue
		else:
			while True:
				random_uuid = uuid.uuid1()
				analysis_id = 'A-%s' % random_uuid.hex[:8]
				db_cur.execute("SELECT analysisID FROM Analysis WHERE analysisID='%s'" % analysis_id)
				if db_cur.fetchone() is None:
					barcodes_json[barcode]['analysis_id'] = analysis_id
					json_modified = True
					break

		db_cur.execute("SELECT panelID FROM Panel WHERE panelID='%s'" % barcodes_json[barcode]['panel'])
		if db_cur.fetchone() is None:
			print "\t*WARNING* (PANEL %s not found in DB)** " % barcodes_json[barcode]['panel']
			print "\t** Panel is needed for foreign key constraint. Update DB first."
			exit()

		db_cur.execute("INSERT INTO Analysis (analysisID, sample, barcode, run, panel, bamPath, analysisDate) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s')" % (analysis_id, barcodes_json[barcode]['sample_id'], barcode, run_name, barcodes_json[barcode]['panel'], barcodes_json[barcode]['bam'], time.strftime("%Y-%m-%d")))
		db_con.commit()

	logging.info("\t - %s, analysisID = %s, panel = %s" % (barcodes_json[barcode]['sample'],barcodes_json[barcode]['analysis_id'],barcodes_json[barcode]['panel']))

if json_modified:
	print "- Re-writing barcodes JSON..."
	json_text = json.dumps(barcodes_json, indent=4, sort_keys=True)
	bc_json = open('%s/barcodes.json' % run_folder,'w')
	bc_json.write(json_text)
	bc_json.close()

###   __   __   ___     __   __   __   __   ___  __   __          __    ###
###  |__) |__) |__  __ |__) |__) /  \ /  ` |__  /__` /__` | |\ | / _`   ###
###  |    |  \ |___    |    |  \ \__/ \__, |___ .__/ .__/ | | \| \__>   ###

if (not options.skip_preprocessing) and (platform != 'ion torrent'):
	logging.info("\n- [%s] PRE-PROCESSING (illumina) ..." % (time.strftime("%H:%M:%S")))
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

###  __   __        ___  __        __   ___                             __     __  ###
### /  ` /  \ \  / |__  |__)  /\  / _` |__      /\  |\ |  /\  |    \ / /__` | /__` ###
### \__, \__/  \/  |___ |  \ /~~\ \__> |___    /~~\ | \| /~~\ |___  |  .__/ | .__/ ###

if not options.skip_coverage_analysis:
	for barcode in ordered_barcodes:
		logging.info("\t\t- [%s] coverageAnalysis ..." % (time.strftime("%H:%M:%S")))
		cmd_list = []

		fastqc_folder = '%s/fastqc' % barcodes_json[barcode]['intermediate_folder']
		coverage_folder = '%s/coverage' % barcodes_json[barcode]['intermediate_folder']
		mosdepth_folder = '%s/mosdepth' % barcodes_json[barcode]['intermediate_folder']

		for f in [fastqc_folder,coverage_folder,mosdepth_folder]:
			if not os.path.isdir(f):
				subprocess.call(['mkdir',f])

		# QC : FASTQC BAM
		logging.info("\t\t- [%s] FastQC (BAM)..." % (time.strftime("%H:%M:%S")))
		cmd = subprocess.Popen(['perl','%s/FastQC/fastqc' % pipeline_folder,'--outdir',fastqc_folder,barcodes_json[barcode]['bam']],stdout=open('%s/fastqc_bam.stdout.txt' % fastqc_folder,'w'),stderr=open('%s/fastqc_bam.stderr.txt' % fastqc_folder,'w'))

		# SAMTOOLS STATS & DEPTH
		logging.info("\t\t- [%s] Samtools stats ..." % (time.strftime("%H:%M:%S")))
		cmd = subprocess.Popen(['samtools', 'stats','-d','-t', barcodes_json[barcode]['target_bed'],barcodes_json[barcode]['bam']],stdout=open('%s/%s_%s.stats.txt' % (barcodes_json[barcode]['intermediate_folder'],barcodes_json[barcode]['sample'],barcode),'w'))
		cmd = subprocess.Popen(['samtools', 'depth','-b',barcodes_json[barcode]['target_bed'],barcodes_json[barcode]['bam']],stdout=open('%s/depth.txt' % coverage_folder,'w'), stderr=open('%s/samtools_depth.stderr.txt' % coverage_folder,'w'))

		# MOSDEPTH
		logging.info("\t\t- [%s] mosdepth ..." % (time.strftime("%H:%M:%S")))
		os.chdir(mosdepth_folder)
		cmd = subprocess.Popen(['%s/mosdepth/mosdepth' % pipeline_folder,'-b', barcodes_json[barcode]['target_bed'],'%s_%s' % (barcodes_json[barcode]['sample'],barcode),barcodes_json[barcode]['bam']],stdout=open('%s/mosdepth.stdout.txt' % mosdepth_folder,'w'), stderr=open('%s/mosdepth.stderr.txt' % mosdepth_folder,'w'))

		# BBCTOOLS
		if platform == 'ion torrent':
			special_param = '-a' # -a = amplicon
		elif platform == 'illumina':
			special_param = '-d' # -d = ignore duplicates (pour target coverage)
		cmd = subprocess.Popen([
			'bash', '%s/coverageAnalysis/run_coverage_analysis.sh' % pipeline_folder,
			'-L', 'hg19',
			special_param,
			'-g',
			'-D',coverage_folder,
			'-B',barcodes_json[barcode]['target_bed'],
			barcodes_json[barcode]['reference'],
			barcodes_json[barcode]['bam']
			],
			stdout=open('%s/run_coverage_analysis.stdout.txt' % coverage_folder,'w'), 
			stderr=open('%s/run_coverage_analysis.stderr.txt' % coverage_folder,'w'))
		cmd.communicate()

#####                ___    __   __  
##### |\/| |  | |     |  | /  \ /  ` 
##### |  | \__/ |___  |  | \__X \__, 

if options.run:
	logging.info(" [%s] MultiQC ..." % time.strftime("%H:%M:%S"))
	multiqc_folder = '%s/_multiqc' % run_folder
	if not os.path.isdir(multiqc_folder):
		subprocess.call(['mkdir', multiqc_folder])
	cmd = subprocess.Popen(['/usr/local/bin/multiqc',
	'-f','-c','%s/multiqc/custom_config.yaml' % pipeline_folder,
	'-o', '%s/_multiqc' % run_folder,
	'-m', 'fastqc', '-m', 'samtools', '-m', 'mosdepth',
	run_folder],
	stdout=open('%s/multiqc.stdout.txt' % multiqc_folder,'w'), 
	stderr=open('%s/multiqc.stderr.txt' % multiqc_folder,'w'))

#####  __        __  ___  __   __        ___  __        __   ___ ###
##### |__) |    /  \  |  /  ` /  \ \  / |__  |__)  /\  / _` |__  ###
##### |    |___ \__/  |  \__, \__/  \/  |___ |  \ /~~\ \__> |___ ###

# if options.run:
	# logging.info(" [%s] plotCoverage ..." % time.strftime("%H:%M:%S"))
	# if not os.path.isdir('%s/_plotCoverage' % run_folder):
		# subprocess.call(['mkdir', '%s/_plotCoverage' % run_folder])
	# cmd = subprocess.Popen(['python', '%s/plotCoverage/plotCoverage.py' % pipeline_folder,'--run-folder', run_folder],stdout=open('%s/_plotCoverage/plotCoverage.stdout.txt' % run_folder,'w'), stderr=open('%s/_plotCoverage/plotCoverage.stderr.txt' % run_folder,'w'))

#####  __            ###
##### /  ` |\ | \  / ###
##### \__, | \|  \/  ###

if options.run:
	logging.info(" [%s] CNV Analysis ..." % (time.strftime("%H:%M:%S")))
	if not os.path.isdir('%s/_CNA' % run_folder):
		subprocess.call(['mkdir', '%s/_CNA' % run_folder])
	cmd = subprocess.Popen(['python','%s/CNV/run_cna.py' % pipeline_folder, run_folder], stdout=open('%s/_CNA/run_cna.stdout.txt' % run_folder,'w'), stderr=open('%s/_CNA/run_cna.stderr.txt' % run_folder,'w'))

#####            __              ___     __                         __                __                     __  ___      ___    __       #####
##### \  /  /\  |__) |  /\  |\ |  |  __ /  `  /\  |    |    | |\ | / _`     /\  |\ | |  \     /\  |\ | |\ | /  \  |   /\   |  | /  \ |\ | #####
#####  \/  /~~\ |  \ | /~~\ | \|  |     \__, /~~\ |___ |___ | | \| \__>    /~~\ | \| |__/    /~~\ | \| | \| \__/  |  /~~\  |  | \__/ | \| #####

logging.info("\n- [%s] VARIANT-CALLING AND ANNOTATION ..." % (time.strftime("%H:%M:%S")))
sampleask = True
for barcode in ordered_barcodes:
	logging.info("\t- %s :" % (barcodes_json[barcode]['sample']))
	if platform == 'illumina':
		mutect2_folder  = '%s/mutect2' % barcodes_json[barcode]['intermediate_folder']
		vardict_folder  = '%s/vardict' % barcodes_json[barcode]['intermediate_folder']
		varscan2_folder = '%s/varscan2' % barcodes_json[barcode]['intermediate_folder']
		lofreq_folder   = '%s/lofreq' % barcodes_json[barcode]['intermediate_folder']

		# FOR SKIPING SAMPLES # WARNING check-contamination will not work #
		# if sampleask:
			# proceed = raw_input('continue? (y/n/stopask)\n')
			# if proceed == 'y' or proceed == 'yes':
				# pass
			# elif proceed == 'stopask':
				# sampleask = False
				# pass
			# else:
				# continue

		if not options.skip_calling:
			# CREATE FOLDERS
			if os.path.isdir(lofreq_folder):
				subprocess.call(['rm','-r',lofreq_folder]) # because "cowardly refusing to overwrite" lofreq bullshit
			for f in [mutect2_folder,vardict_folder,varscan2_folder,lofreq_folder]:
				if not os.path.isdir(f):
					subprocess.call(['mkdir',f])
			# SETTINGS
			mutect2_thread = 12 # x hmm-threads
			vardict_thread = 12
			vardict_ps_list = []
			mutect2_ps_list = []
			vardict_vcf_chunk_list = []
			mutect2_vcf_chunk_list = []

			# VARDICT
			os.chdir(vardict_folder)
			logging.info("\t\t - [%s] VarDict : split_bed.py ..." % time.strftime("%H:%M:%S"))
			cmd = subprocess.Popen(['python','%s/scripts/split_bed.py' % pipeline_folder,'--bed',barcodes_json[barcode]['target_bed'],'--scatter-count',str(vardict_thread),'--output-folder',vardict_folder],stdout=open('%s/split_bed.stdout.txt'%vardict_folder,'w'),stderr=open('%s/split_bed.stderr.txt'%vardict_folder,'w'))
			cmd.communicate()
			logging.info("\t\t - [%s] VarDict : running %s chunks ..." % (time.strftime("%H:%M:%S"),vardict_thread))
			FNULL = open(os.devnull, 'w')
			for i in range(vardict_thread):
				vardict_bed_chunk = '%s/%04d-scattered.bed' % (vardict_folder,i)
				vardict_vcf_chunk = '%s/%s.vardict.%04d.vcf' % (vardict_folder,barcodes_json[barcode]['sample'],i)
				vardict_vcf_chunk_list.append(vardict_vcf_chunk)
				ps = subprocess.Popen(['perl','%s/vardict/VarDictJava/VarDict/vardict' % pipeline_folder,
					'-G', barcodes_json[barcode]['reference'],
					'-f','0.01',
					'-N',barcodes_json[barcode]['sample'],
					'-b', barcodes_json[barcode]['bam'],
					'-c','1','-S','2','-E','3','-g','4',
					vardict_bed_chunk],
					stdout=open('%s/vardict.%04d.stdout.txt' % (vardict_folder,i),'w'),
					stderr=open('%s/vardict.%04d.stderr.txt' % (vardict_folder,i),'w'))
				vardict_ps_list.append(ps)

			# MUTECT2
			logging.info("\t\t - [%s] Mutect2 : SplitIntervals ..." % time.strftime("%H:%M:%S"))
			cmd = subprocess.Popen(['docker','exec','-it','gatk','gatk','SplitIntervals','-R',barcodes_json[barcode]['reference'],'-L',barcodes_json[barcode]['intervals'],'-O',mutect2_folder,'--scatter-count',str(mutect2_thread)],stdout=open('%s/splitintervals.stdout.txt'%mutect2_folder,'w'),stderr=open('%s/splitintervals.stderr.txt'%mutect2_folder,'w'))
			cmd.communicate()
			logging.info("\t\t - [%s] Mutect2 : running %s chunks ..." % (time.strftime("%H:%M:%S"),mutect2_thread))
			FNULL = open(os.devnull, 'w')
			for i in range(mutect2_thread):
				mutect2_interval_chunk = '%s/%04d-scattered.interval_list' % (mutect2_folder,i)
				mutect2_vcf_chunk = '%s/%s.mutect2.chunk%04d.vcf' % (mutect2_folder,barcodes_json[barcode]['sample'],i)
				mutect2_vcf_chunk_list.append(mutect2_vcf_chunk)
				ps = subprocess.Popen(['docker','exec','-it','gatk','gatk','Mutect2',
				'-I',barcodes_json[barcode]['bam'],
				'-R',barcodes_json[barcode]['reference'],
				'-L',mutect2_interval_chunk,
				'--max-mnp-distance','0',
				'--max-reads-per-alignment-start','0',# '--native-pair-hmm-threads','1',
				'-O',mutect2_vcf_chunk],stdout=open('%s/mutect2.chunk%04d.stdout.txt' % (mutect2_folder,i),'w'),stderr=open('%s/mutect2.chunk%04d.stderr.txt' % (mutect2_folder,i),'w'))
				mutect2_ps_list.append(ps)
			subprocess.call(['stty','sane'])

			# VARSCAN2
			logging.info("\t\t - [%s] VarScan2 : mpileup ..." % time.strftime("%H:%M:%S"))
			varscan2_mpileup_ps = subprocess.Popen(['samtools','mpileup','-f',barcodes_json[barcode]['reference'],'-l',barcodes_json[barcode]['target_bed'],'-o','%s/%s.pileup' % (varscan2_folder,barcodes_json[barcode]['sample']),barcodes_json[barcode]['bam']],stdout=open('%s/mpileup.stdout.txt'%varscan2_folder,'w'),stderr=open('%s/mpileup.stderr.txt'%varscan2_folder,'w'))

			# LOFREQ
			logging.info("\t\t - [%s] LoFreq : indelqual ..." % time.strftime("%H:%M:%S"))
			if barcodes_json[barcode]['platform'].lower() == 'illumina':
				lofreq_ps = subprocess.Popen(['%s/lofreq/bin/lofreq' % pipeline_folder,'indelqual','--dindel','-f',barcodes_json[barcode]['reference'],'-o','%s/%s.indelqual.bam' % (lofreq_folder,barcodes_json[barcode]['sample']),barcodes_json[barcode]['bam']], stdout=open('%s/lofreq_indelqual_dindel.stdout.txt' % lofreq_folder,'w'), stderr=open('%s/lofreq_indelqual_dindel.stderr.txt' % lofreq_folder,'w'))
			else: # Set INDELQUAL uniform 45 for ion torrent
				lofreq_ps = subprocess.Popen(['%s/lofreq/bin/lofreq' % pipeline_folder,'indelqual','--uniform','45','-o','%s/%s.indelqual.bam' % (lofreq_folder,barcodes_json[barcode]['sample']),barcodes_json[barcode]['bam']], stdout=open('%s/lofreq_indelqual_uniform.stdout.txt' % lofreq_folder,'w'), stderr=open('%s/lofreq_indelqual_uniform.stderr.txt' % lofreq_folder,'w'))

			calling_in_progress = True
			varscan2_running = True
			lofreq_running = True
			mutect2_running = True
			vardict_running = True
			while calling_in_progress:
				# VARSCAN2
				if varscan2_running and varscan2_mpileup_ps.poll() is not None : # ps.poll() is not None -> process had ended
					logging.info("\t\t - [%s] VarScan2 : mpileup2cns ..." % time.strftime("%H:%M:%S"))
					cmd = subprocess.Popen(['java','-jar','%s/varscan2/VarScan.v2.4.2.jar' % pipeline_folder,'mpileup2cns','%s/%s.pileup' % (varscan2_folder,barcodes_json[barcode]['sample']),'--min-var-freq','0.02','--output-vcf','1','--variants'],stdout=open('%s/%s.varscan2.cns.vcf' % (varscan2_folder,barcodes_json[barcode]['sample']),'w'),stderr=open('%s/varscan2.cns.stderr.txt'%varscan2_folder,'w'))
					cmd.communicate()
					logging.info("\t\t - [%s] VarScan2 : filter ..." % time.strftime("%H:%M:%S"))
					cmd = subprocess.Popen(['java','-jar','%s/varscan2/VarScan.v2.4.2.jar' % pipeline_folder,'filter','%s/%s.varscan2.cns.vcf' % (varscan2_folder,barcodes_json[barcode]['sample']),
						'--min-coverage','10',
						'--min-reads2', '2',
						'--min-strands2','1',
						'--min-avg-qual','10',
						'--min-var-freq','0.01',
						'--output-file','%s/%s.varscan2.filtered.vcf' % (varscan2_folder,barcodes_json[barcode]['sample'])
					],stdout=open('%s/varscan2.filter.stdout.txt'%varscan2_folder,'w'),stderr=open('%s/varscan2.filter.stderr.txt'%varscan2_folder,'w'))
					cmd.communicate()
					varscan2_running = False
				# LOFREQ
				if lofreq_running and lofreq_ps.poll() is not None :
					logging.info("\t\t - [%s] LoFreq : index indelqual ..." % time.strftime("%H:%M:%S"))
					lofreq_ps = subprocess.call(['samtools','index','%s/%s.indelqual.bam' % (lofreq_folder,barcodes_json[barcode]['sample'])], stdout=open('%s/samtools_index.stdout.txt' % lofreq_folder,'w'), stderr=open('%s/samtools_index.stderr.txt' % lofreq_folder,'w'))
					logging.info("\t\t - [%s] LoFreq : call ..." % time.strftime("%H:%M:%S"))
					cmd = subprocess.Popen(['%s/lofreq/bin/lofreq' % pipeline_folder,'call-parallel','--pp-threads','12','--call-indels','-f',barcodes_json[barcode]['reference'],'-l',barcodes_json[barcode]['target_bed'],'-o','%s/%s.lofreq.vcf' % (lofreq_folder,barcodes_json[barcode]['sample']),'%s/%s.indelqual.bam' % (lofreq_folder,barcodes_json[barcode]['sample'])],stdout=open('%s/lofreq.stdout.txt'%lofreq_folder,'w'),stderr=open('%s/lofreq.stderr.txt'%lofreq_folder,'w'))
					cmd.communicate()
					logging.info("\t\t - [%s] LoFreq : filter ..." % time.strftime("%H:%M:%S"))
					cmd = subprocess.Popen(['%s/lofreq/bin/lofreq' % pipeline_folder,'filter',
						'-i','%s/%s.lofreq.vcf' % (lofreq_folder,barcodes_json[barcode]['sample']),
						'-o','%s/%s.lofreq.filtered.vcf' % (lofreq_folder,barcodes_json[barcode]['sample']),
						'--cov-min','10',
						'--af-min','0.01',
						'--snvqual-thresh','10'
					],stdout=open('%s/lofreq.filter.stdout.txt'%lofreq_folder,'w'),stderr=open('%s/lofreq.filter.stderr.txt'%lofreq_folder,'w'))
					cmd.communicate()
					# remove temp indelqual bam
					subprocess.call(['rm','%s/%s.indelqual.bam' % (lofreq_folder,barcodes_json[barcode]['sample'])])
					lofreq_running = False
				# VARDICT
				vardict_ps_states = [ps.poll() is not None for ps in vardict_ps_list]
				if vardict_running and (False not in vardict_ps_states) :
					logging.info("\t\t - [%s] VarDict : teststrandbias.R ..." % time.strftime("%H:%M:%S"))
					ps = []
					for i in range(vardict_thread):
						cmd = subprocess.Popen(['%s/vardict/VarDictJava/VarDict/teststrandbias.R' % pipeline_folder],
						stdin=open('%s/vardict.%04d.stdout.txt' % (vardict_folder,i),'r'),
						stdout=open('%s/teststrandbias.%04d.stdout.txt' % (vardict_folder,i),'w'),
						stderr=open('%s/teststrandbias.%04d.stderr.txt' % (vardict_folder,i),'w'))
						ps.append(cmd)
					for cmd in ps:
						cmd.wait()
					logging.info("\t\t - [%s] VarDict : var2vcf_valid.pl ..." % time.strftime("%H:%M:%S"))
					ps = []
					for i in range(vardict_thread):
						cmd = subprocess.Popen(['perl','%s/vardict/VarDictJava/VarDict/var2vcf_valid.pl' % pipeline_folder,
							'-N',barcodes_json[barcode]['sample'],
							'-E', # If set, do not print END tag
							'-A', # output all variants at same position
							'-f','0.01'],
						stdin=open('%s/teststrandbias.%04d.stdout.txt' % (vardict_folder,i),'r'),
						stdout=open('%s/%s.vardict.%04d.vcf' % (vardict_folder,barcodes_json[barcode]['sample'],i),'w'),
						stderr=open('%s/var2vcf_valid.%04d.stderr.txt' % (vardict_folder,i),'w'))
						ps.append(cmd)
					for cmd in ps:
						cmd.wait()
					logging.info("\t\t - [%s] VarDict : UpdateVCFSequenceDictionary ..." % time.strftime("%H:%M:%S"))
					ps = []
					for vardict_vcf_chunk in vardict_vcf_chunk_list:
						cmd = subprocess.Popen(['docker','exec','-it','gatk','gatk','UpdateVCFSequenceDictionary',
							'-V',vardict_vcf_chunk,
							'--source-dictionary',barcodes_json[barcode]['bam'],
							'--output',vardict_vcf_chunk.replace('.vcf','.contig.vcf')],
						stdout=open('%s/UpdateVCFSequenceDictionary.stdout.txt' % vardict_folder,'w'),
						stderr=open('%s/UpdateVCFSequenceDictionary.stderr.txt' % vardict_folder,'w'))
						ps.append(cmd)
					for cmd in ps:
						cmd.wait()
					subprocess.call(['stty','sane'])
					logging.info("\t\t - [%s] VarDict : GatherVcfs ..." % time.strftime("%H:%M:%S"))
					args = ['docker','exec','-it','gatk','gatk','GatherVcfs','-O','%s/%s.vardict.vcf' % (vardict_folder,barcodes_json[barcode]['sample'])]
					for vardict_vcf_chunk in vardict_vcf_chunk_list:
						args.append('-I')
						args.append(vardict_vcf_chunk.replace('.vcf','.contig.vcf'))
					cmd = subprocess.Popen(args,stdout=open('%s/gathervcf.stdout.txt' % vardict_folder,'w'),stderr=open('%s/gathervcf.stderr.txt' % vardict_folder,'w'))
					cmd.wait()
					subprocess.call(['stty','sane'])
					vardict_running = False
				# MUTECT2
				mutect2_ps_states = [ps.poll() is not None for ps in mutect2_ps_list]
				if mutect2_running and (False not in mutect2_ps_states) :
					logging.info("\t\t - [%s] Mutect2 : GatherVcfs ..." % time.strftime("%H:%M:%S"))
					args = ['docker','exec','-it','gatk','gatk','GatherVcfs','-O','%s/%s.mutect2.vcf' % (mutect2_folder,barcodes_json[barcode]['sample'])]
					for mutect2_vcf_chunk in mutect2_vcf_chunk_list:
						args.append('-I')
						args.append(mutect2_vcf_chunk)
					cmd = subprocess.Popen(args,stdout=open('%s/gathervcf.stdout.txt' % mutect2_folder,'w'),stderr=open('%s/gathervcf.stderr.txt' % mutect2_folder,'w'))
					cmd.wait()
					subprocess.call(['stty','sane'])
					logging.info("\t\t - [%s] Mutect2 : MergeMutectStats ..." % time.strftime("%H:%M:%S"))
					args = ['docker','exec','-it','gatk','gatk','MergeMutectStats','-O','%s/%s.mutect2.vcf.stats' % (mutect2_folder,barcodes_json[barcode]['sample'])]
					for mutect2_vcf_chunk in mutect2_vcf_chunk_list:
						args.append('-stats')
						args.append(mutect2_vcf_chunk+'.stats')
					cmd = subprocess.Popen(args,stdout=open('%s/mergemutectstats.stdout.txt' % mutect2_folder,'w'),stderr=open('%s/mergemutectstats.stderr.txt' % mutect2_folder,'w'))
					cmd.wait()
					subprocess.call(['stty','sane'])
					logging.info("\t\t - [%s] Mutect2 : FilterMutectCalls ..." % time.strftime("%H:%M:%S"))
					cmd = subprocess.Popen([
						'docker','exec','-it','gatk','gatk', 'FilterMutectCalls',
						'-V','%s/%s.mutect2.vcf' % (mutect2_folder,barcodes_json[barcode]['sample']),
						'-R',barcodes_json[barcode]['reference'],
						'--max-events-in-region', '20',
						'--min-allele-fraction', '0.01',
						'--min-median-base-quality','10',
						'--unique-alt-read-count','10',
						'-O','%s/%s.mutect2.filtered.vcf' % (mutect2_folder,barcodes_json[barcode]['sample'])
						#'--panel-of-normals', '%s/reference_files/mutect/colon_lung_pon.vcf.gz' % pipeline_folder
					], stdout=open('%s/filtermutectcalls.stdout.txt' % mutect2_folder,'w'),stderr=open('%s/filtermutectcalls.stderr.txt' % mutect2_folder,'w'))
					cmd.communicate()
					subprocess.call(['stty','sane'])
					mutect2_running = False

				subprocess.call(['stty','sane'])
				if not varscan2_running and not lofreq_running and not vardict_running and not mutect2_running:
					calling_in_progress = False
				else:
					time.sleep(5)
	elif platform == 'ion torrent':
		tvc_de_novo_folder = '%s/tvc_de_novo' % barcodes_json[barcode]['intermediate_folder']
		tvc_only_hotspot_folder = '%s/tvc_only_hotspot' % barcodes_json[barcode]['intermediate_folder']
		subprocess.call(['mkdir',tvc_de_novo_folder])
		subprocess.call(['mkdir',tvc_only_hotspot_folder])

		# RUN 1 : DE NOVO
		logging.info("\t - [%s] variantCaller <de novo> ..." % (time.strftime("%H:%M:%S")))
		cmd = subprocess.Popen([
			'python', '%s/variantCaller/bin/variant_caller_pipeline.py' % pipeline_folder,
			'--input-bam',       barcodes_json[barcode]['bam'],
			'--output-dir',      '%s/tvc_de_novo' % barcodes_json[barcode]['intermediate_folder'],
			'--reference-fasta', barcodes_json[barcode]['reference'],
			'--region-bed',      barcodes_json[barcode]['merged_bed'],
			'--primer-trim-bed', barcodes_json[barcode]['target_bed'],
			'--parameters-file', barcodes_json[barcode]['tvc_param'],
			'--error-motifs', '%s/variantCaller/share/TVC/sse/motifset.txt' % pipeline_folder,
			'--postprocessed-bam','%s/tvc_de_novo/processed.bam' % barcodes_json[barcode]['intermediate_folder'],
			], 
			stdout=open('%s/tvc_de_novo/variant_caller_pipeline.stdout.txt' % barcodes_json[barcode]['intermediate_folder'],'w'),
			stderr=open('%s/tvc_de_novo/variant_caller_pipeline.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'w'))
		cmd.communicate()

		subprocess.call(['gzip','-d','-c','%s/tvc_de_novo/TSVC_variants.vcf.gz' % barcodes_json[barcode]['intermediate_folder']],stdout=open('%s/tvc_de_novo/TSVC_variants.vcf' % barcodes_json[barcode]['intermediate_folder'],'w'))

		# Generate Variant Table ('alleles.xls' file)
		cmd = subprocess.Popen([
			'python', '%s/variantCaller/bin/generate_variant_tables.py' % pipeline_folder,
			'--input-vcf',		'%s/tvc_de_novo/TSVC_variants.vcf' % barcodes_json[barcode]['intermediate_folder'],
			'--region-bed',		barcodes_json[barcode]['target_bed'],
			'--output-xls',		'%s/tvc_de_novo/output.xls' % barcodes_json[barcode]['intermediate_folder'],
			'--alleles2-xls',	'%s/tvc_de_novo/alleles.xls' % barcodes_json[barcode]['intermediate_folder']
			], 
			stdout=open('%s/tvc_de_novo/generate_variant_tables.stdout.txt' % barcodes_json[barcode]['intermediate_folder'],'w'), 
			stderr=open('%s/tvc_de_novo/generate_variant_tables.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'w'))
		cmd.communicate()
		
		# RUN 2 : ONLY HOTSPOT
		logging.info("\t - [%s] variantCaller <only hotspot> ..." % (time.strftime("%H:%M:%S")))
		if barcodes_json[barcode]['hotspot_vcf'] != '':
			cmd = subprocess.Popen([
				'python', '%s/variantCaller/bin/variant_caller_pipeline.py' % pipeline_folder,
				'--input-bam',       barcodes_json[barcode]['bam'],
				'--output-dir',      '%s/tvc_only_hotspot' % barcodes_json[barcode]['intermediate_folder'],
				'--reference-fasta', barcodes_json[barcode]['reference'],
				'--region-bed',      barcodes_json[barcode]['merged_bed'],
				'--primer-trim-bed', barcodes_json[barcode]['target_bed'],
				'--parameters-file', barcodes_json[barcode]['tvc_param_hotspot'],
				'--error-motifs',    '%s/variantCaller/share/TVC/sse/motifset.txt' % pipeline_folder,
				'--hotspot-vcf',     barcodes_json[barcode]['hotspot_vcf'],
				], 
				stdout=open('%s/tvc_only_hotspot/variant_caller_pipeline.stdout.txt' % barcodes_json[barcode]['intermediate_folder'],'w'),
				stderr=open('%s/tvc_only_hotspot/variant_caller_pipeline.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'w'))
			cmd.communicate()

			subprocess.call(['gzip','-d','-c','%s/tvc_only_hotspot/TSVC_variants.vcf.gz' % barcodes_json[barcode]['intermediate_folder']],stdout=open('%s/tvc_only_hotspot/TSVC_variants.vcf' % barcodes_json[barcode]['intermediate_folder'],'w'))

			# Generate Variant Table ('alleles.xls' file)
			cmd = subprocess.Popen([
				'python', '%s/variantCaller/bin/generate_variant_tables.py' % pipeline_folder,
				'--input-vcf',		'%s/tvc_only_hotspot/TSVC_variants.vcf' % barcodes_json[barcode]['intermediate_folder'],
				'--region-bed',		barcodes_json[barcode]['target_bed'],
				'--hotspots',
				'--output-xls',		'%s/tvc_only_hotspot/output.xls' % barcodes_json[barcode]['intermediate_folder'],
				'--alleles2-xls',	'%s/tvc_only_hotspot/alleles.xls' % barcodes_json[barcode]['intermediate_folder']
				], 
				stdout=open('%s/tvc_only_hotspot/generate_variant_tables.stdout.txt' % barcodes_json[barcode]['intermediate_folder'],'w'),
				stderr=open('%s/tvc_only_hotspot/generate_variant_tables.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'w'))
			cmd.communicate()

		###         __   ___  __  ___     __   __  ### 
		### | |\ | /__` |__  |__)  |     |  \ |__) ### 
		### | | \| .__/ |___ |  \  |     |__/ |__) ### 

	if not options.skip_insert_db:
		logging.info("\t\t - [%s] insert variants and metrics into DB ..." % (time.strftime("%H:%M:%S")))
		abl1 = 'no'
		if 'ABL1_NM_005157.fasta' in barcodes_json[barcode]['reference']:
			abl1 = 'yes'

		cmd = subprocess.Popen(['python','%s/variantBase/insert_db_variants_illumina.py' % pipeline_folder,
			'--analysis', barcodes_json[barcode]['analysis_id'],
			'--vcf', '%s/%s.mutect2.filtered.vcf' % (mutect2_folder,barcodes_json[barcode]['sample']),
			'--vcf', '%s/%s.varscan2.filtered.vcf' % (varscan2_folder,barcodes_json[barcode]['sample']),
			'--vcf', '%s/%s.lofreq.filtered.vcf' % (lofreq_folder,barcodes_json[barcode]['sample']),
			'--vcf', '%s/%s.vardict.vcf' % (vardict_folder,barcodes_json[barcode]['sample']),
			'--abl1', abl1
			],
			stdout=open('%s/insert_db_variants.stdout.txt' % barcodes_json[barcode]['intermediate_folder'],'w'),
			stderr=open('%s/insert_db_variants.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'w'))
		cmd.communicate()
		with open('%s/insert_db_variants.stdout.txt' % barcodes_json[barcode]['intermediate_folder'],'r') as stdout:
			for line in stdout.readlines():
				print '\t\t\t' + line.replace('\n','')
		with open('%s/insert_db_variants.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'r') as stderr:
			for line in stderr.readlines():
				print '\t\t\t' + line.replace('\n','')


	if not options.skip_annotation:
		###                 __  ___      ___  ___          ___                    __              ___  __  ###
		###  /\  |\ | |\ | /  \  |   /\   |  |__     |\ | |__  |  |    \  /  /\  |__) |  /\  |\ |  |  /__` ###
		### /~~\ | \| | \| \__/  |  /~~\  |  |___    | \| |___ |/\|     \/  /~~\ |  \ | /~~\ | \|  |  .__/ ###

		logging.info("\t\t - [%s] annotate new DB variants ..." % (time.strftime("%H:%M:%S")))

		logging.info("\t\t\t - [%s] HGVS check ..." % (time.strftime("%H:%M:%S")))
		cmd = subprocess.Popen(['python','%s/variantAnnotation/annotate_variantbase.step0.py' % pipeline_folder,'--new'], #--new only annote new variants
		stdout=open('%s/annotate_variantbase.step0.stdout.txt' % barcodes_json[barcode]['intermediate_folder'],'w'), 
		stderr=open('%s/annotate_variantbase.step0.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'w'))
		cmd.communicate()
		with open('%s/annotate_variantbase.step0.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'r') as stderr:
			for line in stderr.readlines():
				if not line.startswith('No handlers could be found for logger "hgvs"'):
					print '\t\t\t' + line.replace('\n','')

		logging.info("\t\t\t - [%s] Annovar and VEP ..." % (time.strftime("%H:%M:%S")))
		cmd = subprocess.Popen(['python','%s/variantAnnotation/annotate_variantbase.step1.py' % pipeline_folder,'--new','--output-folder',barcodes_json[barcode]['intermediate_folder']], 
		stdout=open('%s/annotate_variantbase.step1.stdout.txt' % barcodes_json[barcode]['intermediate_folder'],'w'), 
		stderr=open('%s/annotate_variantbase.step1.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'w'))
		cmd.communicate()
		with open('%s/annotate_variantbase.step1.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'r') as stderr:
			for line in stderr.readlines():
				print '\t\t\t' + line.replace('\n','')

		logging.info("\t\t\t - [%s] Merging annotations and updating DB ..." % (time.strftime("%H:%M:%S")))
		if os.path.isfile('%s/annovar/annovar_input.tsv' % barcodes_json[barcode]['intermediate_folder']):
			cmd = subprocess.Popen(['python','%s/variantAnnotation/annotate_variantbase.step2.py' % pipeline_folder,'--annovar-results', '%s/annovar/annovar.hg19_multianno.txt' % barcodes_json[barcode]['intermediate_folder'],'--vep-results', '%s/vep/vep_output.tsv' % barcodes_json[barcode]['intermediate_folder']], 
			stdout=open('%s/annotate_variantbase.step2.stdout.txt' % barcodes_json[barcode]['intermediate_folder'],'w'),
			stderr=open('%s/annotate_variantbase.step2.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'w'))
			cmd.communicate()
			with open('%s/annotate_variantbase.step2.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'r') as stderr:
				for line in stderr.readlines():
					print '\t\t\t' + line.replace('\n','')

###  ___                   __   ___  __   __   __  ___ ###
### |__  | |\ |  /\  |    |__) |__  |__) /  \ |__)  |  ###
### |    | | \| /~~\ |___ |  \ |___ |    \__/ |  \  |  ###

	if not os.path.isdir('%s/finalReport' % barcodes_json[barcode]['intermediate_folder']):
		subprocess.call(['mkdir', '%s/finalReport' % barcodes_json[barcode]['intermediate_folder']])
	barcodes_json[barcode]['finalReport'] = '%s/%s_%s.finalReport.xlsx' % (barcodes_json[barcode]['sample_folder'],barcodes_json[barcode]['sample'],barcode)

	xmin = '300' # for coverage analysis sheet in finalreport
	# if 'cDNA' in param_hotspot_only:
		# xmin = '2000'

	logging.info("\t\t - [%s] finalReport ..." % (time.strftime("%H:%M:%S")))
	cmd = subprocess.Popen([
		'python','%s/finalReport/finalReport.py' % pipeline_folder,
		'--analysis', barcodes_json[barcode]['analysis_id'],
		'--xmin', xmin,
	], 
	stdout=open('%s/finalReport/finalReport.stdout.txt' % barcodes_json[barcode]['intermediate_folder'],'w'), 
	stderr=open('%s/finalReport/finalReport.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'w'))
	cmd.communicate()
	with open('%s/finalReport/finalReport.stdout.txt' % barcodes_json[barcode]['intermediate_folder'],'r') as stdout:
		for line in stdout.readlines():
			print '\t\t\t' + line.replace('\n','')
	with open('%s/finalReport/finalReport.stderr.txt' % barcodes_json[barcode]['intermediate_folder'],'r') as stderr:
		for line in stderr.readlines():
			print '\t\t\t' + line.replace('\n','')

if options.run:
	###  __   __   __     __  ___  __  ###
	### /__` /  ` |__) | |__)  |  /__` ###
	### .__/ \__, |  \ | |     |  .__/ ###

	for barcode in ordered_barcodes:
		if 'ACROMETRIX' in barcodes_json[barcode]['sample'].upper():
			subprocess.call(['python','%s/scripts/Acrometrix_check.py' % pipeline_folder,barcodes_json[barcode]['finalReport'],barcodes_json[barcode]['sample'],run_name])
		elif 'HORIZON' in barcodes_json[barcode]['sample'].upper():
			subprocess.call(['python','%s/scripts/Horizon_check.py' % pipeline_folder,barcodes_json[barcode]['finalReport'],barcodes_json[barcode]['sample'],run_name])
		elif 'BARBI' in barcodes_json[barcode]['sample'].upper():
			subprocess.call(['python','%s/scripts/Temoin_TP53_check.py' % pipeline_folder,barcodes_json[barcode]['finalReport'],barcodes_json[barcode]['sample'],run_name])
		for c in ['BAF-','BAF5','P190']:
			if c in barcodes_json[barcode]['sample'].upper():
				subprocess.call(['python','%s/scripts/Temoins_ABL1_check.py' % pipeline_folder,barcodes_json[barcode]['finalReport'],barcodes_json[barcode]['sample'],run_name])
		# if run_type == 'SBT':
			# subprocess.call(['python','%s/scripts/collect_variants_MET_intron_13-14.py' % pipeline_folder,barcodes_json[barcode]['finalReport'],sample,run_name])

	if os.path.isfile('%s/tvc_de_novo/processed.bam' % barcodes_json[barcode]['intermediate_folder']):
		subprocess.call(['mv','%s/tvc_de_novo/processed.bam' % barcodes_json[barcode]['intermediate_folder'],'%s/%s_%s.processed.bam' % (barcodes_json[barcode]['sample_folder'],barcodes_json[barcode]['sample'],barcode)])
	if os.path.isfile('%s/tvc_de_novo/processed.bam.bai' % barcodes_json[barcode]['intermediate_folder']):
		subprocess.call(['mv','%s/tvc_de_novo/processed.bam.bai' % barcodes_json[barcode]['intermediate_folder'],'%s/%s_%s.processed.bam.bai' % (barcodes_json[barcode]['sample_folder'],barcodes_json[barcode]['sample'],barcode)])

	# shutil.make_archive(intermediate_folder,'zip',intermediate_folder)
	# shutil.rmtree(intermediate_folder)

	###  __        ___  __           __   __       ___                        ___    __       ###
	### /  ` |__| |__  /  ` |__/    /  ` /  \ |\ |  |   /\   |\/| | |\ |  /\   |  | /  \ |\ | ###
	### \__, |  | |___ \__, |  \    \__, \__/ | \|  |  /~~\  |  | | | \| /~~\  |  | \__/ | \| ###

	# checkconta_folder = '%s/_checkContamination' % run_folder
	# if not os.path.isdir(checkconta_folder):
		# subprocess.call(['mkdir', checkconta_folder])

	# logging.info(" [%s] checkContamination ..." % time.strftime("%H:%M:%S"))
	# cmd = subprocess.Popen(['python','%s/checkContamination/checkContamination.py' % pipeline_folder,'--run-folder', run_folder],stdout=open('%s/checkContamination.stdout.txt' % checkconta_folder,'w'),stderr=open('%s/checkContamination.stderr.txt' % checkconta_folder,'w'))
	# cmd.communicate()
	# check_stderr('%s/checkContamination.stderr.txt' % checkconta_folder,indent=1)

	####  __        ___  __                  ___ ###
	#### /  ` |__| |__  /  ` |__/  |\/| |  |  |  ###
	#### \__, |  | |___ \__, |  \  |  | \__/  |  ###

	checkMut_folder = '%s/_checkMut' % run_folder
	if not os.path.isdir(checkMut_folder):
		subprocess.call(['mkdir',checkMut_folder])
	logging.info(" [%s] checkMut (routine variants) ..." % time.strftime("%H:%M:%S"))
	subprocess.call(['python','%s/checkMut/routine_checkMut_illumina.py' % pipeline_folder,run_folder])

	#  ___      __   ___                    __              ___  __        __   ___ 
	# |__  \_/ /  ` |__  |       \  /  /\  |__) |  /\  |\ |  |  |__)  /\  /__` |__  
	# |___ / \ \__, |___ |___     \/  /~~\ |  \ | /~~\ | \|  |  |__) /~~\ .__/ |___ 

	subprocess.call(['python', '%s/variantBase/make_excel_VariantBase.py' % pipeline_folder, '--panels', 'LAM-illumina-v1,LAM-illumina-v2', '--output', '/media/n06lbth/sauvegardes_pgm/LAM/panel-capture/VariantBase_capture.xlsx'])
