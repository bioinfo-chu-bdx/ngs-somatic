#!/usr/bin/python
import os
import sys
import time
import glob
import json
import uuid
import shutil
import zipfile
import sqlite3
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


### GATHERING PARAMETERS ############################################################

FNULL = open(os.devnull, 'w')
parser = OptionParser()
parser.add_option('-f', '--fastq',	help="input fastq R1 (auto detect R2 if present)", 	dest='fastq')
parser.add_option('-r', '--run',	help="Run folder path  for FULL RUN ANALYSIS",		dest='run') 
(options, args) = parser.parse_args()

if options.fastq and options.run:
	sys.stderr.write("[run_analysis.py] Error: <--fastq> and <--full-run> are not compatibles\n")
	sys.exit()
if options.run:
	fastqlist = glob.glob('%s/*/*_R1_001.fastq*' % options.run) # .gz?
	run_folder = options.run
elif options.fastq:
	fastqlist = [options.fastq]
	run_folder = os.path.dirname(os.path.dirname(options.fastq))
else:
	sys.stderr.write("[run_analysis.py] Error: no <--fastq> or <--full-run> specified\n")
	sys.exit()

# set up logging to file
if options.run:
	logging.basicConfig(level=logging.INFO,format='%(levelname)-7s %(message)s',filename='/%s/run_analysis.log' % run_folder,filemode='w')
else:
	logging.basicConfig(level=logging.INFO,format='%(levelname)-7s %(message)s',filename='/%s/sample_analysis.log' % os.path.dirname(options.fastq),filemode='w')
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console) # add the handler to the root logger

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

if os.path.isfile('%s/barcodes.json' % run_folder):
	with open('%s/barcodes.json' % run_folder, 'r') as g:
		barcodes_json = json.load(g)
else:
	logging.error('barcodes.json file not found')
	sys.exit() 

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

control_names = ['H2O','H20','NTC'] # liste des noms possibles pour les temoins negatifs
sampleask = False
checkconta_bamlist = []
fastq_data = {}

if run_folder.endswith('/'):
	run_name = os.path.basename(os.path.dirname(run_folder))
else:
	run_name = os.path.basename(run_folder)

######################################################################################

logging.info(" [%s] STARTING ANALYSIS ..." % (time.strftime("%H:%M:%S")))

# CREATE RUN DB ENTRY - IF INEXISTANT -
if options.run:
	db_cur.execute("SELECT runID FROM Run WHERE runID='%s'"%run_name)
	if db_cur.fetchone() is None:
		logging.info("\t - [%s] Create new Run DB entry ..." % (time.strftime("%H:%M:%S")))
		platform = 'NextSeq'
		run_date = date.today()
		try:
			db_cur.execute("INSERT INTO Run (runID, platform, runPath, runDate) VALUES ('%s', '%s', '%s', '%s')" % (run_name,platform,run_folder,run_date))
			db_con.commit()
		except Exception as e:
			logging.warning("*WARNING* (RUN table)** %s" % e)
			sys.exit()

### RUN ANALYSIS ######################################################################

for fastqfile in sorted(fastqlist) :
	sample = fastqfile.split('/')[-1].split('_S')[0]
	barcode = 'S%s' % fastqfile.split('_S')[-1].split('_R')[0]
	sample_folder = os.path.dirname(fastqfile)

	## CREATE SAMPLE DB ENTRY - IF INEXISTANT -
	sample_id = barcodes_json[barcode]['sample_id']
	name = sample.replace(sample_id,'')
	while name.endswith('-') or name.endswith('_'):
		name = name[:-1]
	db_cur.execute("SELECT sampleID FROM Sample WHERE sampleID='%s'" % sample_id)
	if db_cur.fetchone() is None:
		logging.info("\t - [%s] new Sample and Analysis DB entry for %s ..." % (time.strftime("%H:%M:%S"), sample))
		try:
			db_cur.execute("INSERT INTO Sample (sampleID, sampleName) VALUES ('%s', '%s')" % (sample_id,name))
			db_con.commit()
		except Exception as e:
			logging.warning("\t*WARNING* (SAMPLE table)** %s" % e)

	run_type = None
	bed_name = barcodes_json[barcode]['target_region_filepath'].split('/unmerged/detail/')[-1]
	for rt in global_param['run_type']:
		if global_param['run_type'][rt]['target_bed'].split('/')[-1] == bed_name:
			run_type = rt
			break

	## CREATE ANALYSIS DB ENTRY - IF INEXISTANT -
	db_cur.execute("SELECT analysisID FROM Analysis WHERE sample='%s' and run='%s' and panel='%s'"% (sample_id,run_name,bed_name))
	db_analysis = db_cur.fetchone()
	if db_analysis is None:
		try:
			db_cur.execute("SELECT panelID FROM Panel WHERE panelID='%s'" % bed_name)
			if db_cur.fetchone() is None:
				logging.warning("\t*WARNING* (PANEL %s not found in DB)** " % panel)
				logging.info("\t** Panel is needed for foreign key constraint. Update DB first.")
				sys.exit()
			random_uuid = uuid.uuid1()
			analysis_id = 'A-%s' % random_uuid.hex[:8]
			db_cur.execute("INSERT INTO Analysis (analysisID, sample, barcode, run, panel, bamPath, analysisDate) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s')" % (analysis_id, sample_id, barcode, run_name, bed_name, fastqfile, time.strftime("%Y-%m-%d")))
			db_con.commit()
		except Exception as e:
			logging.warning("\t*WARNING* (Analysis table)** %s" % e)
	else:
		analysis_id = db_analysis['analysisID']
	logging.info("\t - %s, analysisID = %s, panel = %s" % (sample,analysis_id,run_type))

	### RUN TYPE PARAMETERS ###
	if not run_type:
		logging.warning("\t -- Warning : run type (panel) not found for %s. Sample will not be processed." % sample)
		continue
	reference = global_param['run_type'][run_type]['reference']
	covered_bed = global_param['run_type'][run_type].get('covered_bed','target_bed') # si covered existe, sinon utiliser le target
	target_bed = global_param['run_type'][run_type]['target_bed']
	if covered_bed == target_bed:
		logging.info("WARNING : covered_bed is also roi_bed")
	merged_bed = global_param['run_type'][run_type]['merged_bed']
	intervals = global_param['run_type'][run_type].get('intervals',False)
	if not intervals:
		logging.info("WARNING intervals file missing, cannot analyse this fastq")
		continue
	param = global_param['run_type'][run_type]['vc_parameters']
	param_hotspot_only = global_param['run_type'][run_type]['vc_parameters_hotspot_only']
	if 'PL' in sample_id : # cDNA sample
		if 'vc_parameters_hotspot_cdna' in global_param['run_type'][run_type]:
			logging.info("- Sample %s is cDNA" % sample)
			param_hotspot_only = global_param['run_type'][run_type]['vc_parameters_hotspot_cdna']
	hotspot_vcf = global_param['run_type'][run_type]['hotspot_vcf']

	# DOSSIER INTERMEDIAIRE
	intermediate_folder = '%s/intermediate_files' % sample_folder
	if os.path.isfile('%s.zip' % intermediate_folder): # UNZIP ALL intermediates_files.zip (if RERUN)
		fnull = open(os.devnull, 'w')
		subprocess.call(['unzip', '-o', '%s.zip' % intermediate_folder,'-d', intermediate_folder],stdout=fnull)
	elif not os.path.isdir(intermediate_folder):
		subprocess.call(['mkdir', intermediate_folder])

	fastq_data[fastqfile] = {
	'barcode':barcode,'sample':sample,'sample_id':sample_id,'analysis_id':analysis_id,'run_type':run_type,'reference':reference,
	'sample_folder':sample_folder,'intermediate_folder':intermediate_folder,'target_bed':target_bed,'intervals':intervals,'merged_bed':merged_bed,
	'covered_bed':covered_bed,'param':param,'param_hotspot_only':param_hotspot_only,'hotspot_vcf':hotspot_vcf}

db_con.close()

###   __   __   ___     __   __   __   __   ___  __   __          __    ###
###  |__) |__) |__  __ |__) |__) /  \ /  ` |__  /__` /__` | |\ | / _`   ###
###  |    |  \ |___    |    |  \ \__/ \__, |___ .__/ .__/ | | \| \__>   ###

logging.info("\n- [%s] PRE-PROCESSING ..." % (time.strftime("%H:%M:%S")))

for fastqfile in fastq_data:
	logging.info("\t- %s :" % (fastq_data[fastqfile]['sample']))
	prep_folder = '%s/pre-processing' % fastq_data[fastqfile]['intermediate_folder']
	if not os.path.isdir(prep_folder):
		subprocess.call(['mkdir', prep_folder])

	fastq_r1 = fastqfile
	paired_end_classic = False
	paired_end_mbc = False
	# RUN PAIRED END WITH MBC
	if os.path.isfile(fastqfile.replace('R1_001','R3_001')):
		paired_end_mbc = True
		fastq_data[fastqfile]['fastq_r2'] = fastqfile.replace('R1_001','R2_001')
		fastq_data[fastqfile]['fastq_r3'] = fastqfile.replace('R1_001','R3_001')
	# RUN PAIRED END NO MBC
	elif os.path.isfile(fastqfile.replace('R1_001','R2_001')):
		paired_end_classic = True
		fastq_data[fastqfile]['fastq_r2'] = fastqfile.replace('R1_001','R2_001')

	# (OPTIONAL) AGENT TRIM
	if paired_end_mbc:
		logging.info("\t\t- [%s] Trim MBC ..." % time.strftime("%H:%M:%S"))
		# cmd = subprocess.Popen(['%s/agent/agent.sh' % pipeline_folder,'trim','-xt','-fq1',fastq_r1,'-fq2',fastq_data[fastqfile]['fastq_r3'],'-out_loc',prep_folder], stdout=open('%s/trim.stdout.txt' % prep_folder,'w'), stderr=open('%s/trim.stderr.txt' % prep_folder,'w'))
		# cmd.communicate()
		fastq_r1 = 'zobzob'
		fastq_data[fastqfile]['fastq_r3'] = 'zibzib'

	# # QC : FASTQC
	# logging.info("\t\t- [%s] FastQC ..." % (time.strftime("%H:%M:%S")))
	# fastqc_folder = '%s/fastqc' % fastq_data[fastqfile]['intermediate_folder']
	# if not os.path.isdir(fastqc_folder):
		# subprocess.call(['mkdir', fastqc_folder])
	# cmd = subprocess.Popen(['perl','%s/FastQC/fastqc' % fastq_r1,'--outdir',fastqc_folder],stdout=open('%s/fastqc_r1.stdout.txt' % fastqc_folder,'w'),stderr=open('%s/fastqc_r1.stderr.txt' % fastqc_folder,'w'))
	# cmd = subprocess.Popen(['perl','%s/FastQC/fastqc' % fastq_data[fastqfile]['fastq_r3'],'--outdir',fastqc_folder],stdout=open('%s/fastqc_r3.stdout.txt' % fastqc_folder,'w'),stderr=open('%s/fastqc_r3.stderr.txt' % fastqc_folder,'w'))

	# ALIGNMENT - FASTQ TO SAM
	logging.info("\t\t- [%s] BWA-MEM Alignment ..." % time.strftime("%H:%M:%S"))
	read_group = '@RG\\tID:%s\\tLB:%s\\tPL:%s\\tPU:%s\\tSM:%s' % ('%s' % fastq_data[fastqfile]['sample_id'],'Target-Myeloid-v1','illumina','NDX550372',fastq_data[fastqfile]['sample']) # READ GROUP (which reflects which library a read belongs to and what lane it was sequenced in on the flowcell)
	# if paired_end_mbc :
		# cmd = subprocess.Popen(['%s/bwa/bwa' % pipeline_folder,'mem','-C','-t','12','-M','-R',read_group,'-v','1','-o','%s/%s.sam' % (prep_folder,fastq_data[fastqfile]['sample']),fastq_data[fastqfile]['reference'],fastq_r1,fastq_data[fastqfile]['fastq_r3']], stdout=open('%s/bwa_mem.stdout.txt' % prep_folder,'w'), stderr=open('%s/bwa_mem.stderr.txt' % prep_folder,'w'))
	# elif paired_end_classic:
		# cmd = subprocess.Popen(['%s/bwa/bwa' % pipeline_folder,'mem','-C','-t','12','-M','-R',read_group,'-v','1','-o','%s/%s.sam' % (prep_folder,fastq_data[fastqfile]['sample']),fastq_data[fastqfile]['reference'],fastq_r1,fastq_data[fastqfile]['fastq_r2']], stdout=open('%s/bwa_mem.stdout.txt' % prep_folder,'w'), stderr=open('%s/bwa_mem.stderr.txt' % prep_folder,'w'))
	# else:
		# cmd = subprocess.Popen(['%s/bwa/bwa' % pipeline_folder,'mem','-C','-t','12','-M','-R',read_group,'-v','1','-o','%s/%s.sam' % (prep_folder,fastq_data[fastqfile]['sample']),fastq_data[fastqfile]['reference'],fastq_r1], stdout=open('%s/bwa_mem.stdout.txt' % prep_folder,'w'), stderr=open('%s/bwa_mem.stderr.txt' % prep_folder,'w'))
	# cmd.communicate()

	bam = '%s/%s_%s.bam' % (fastq_data[fastqfile]['sample_folder'],fastq_data[fastqfile]['sample'],fastq_data[fastqfile]['barcode'])
	if paired_end_classic :
		## MARKDUPLICATESSPARK (& SORT)
		logging.info("\t\t- [%s] Marking duplicates and sorting BAM ..." % time.strftime("%H:%M:%S"))
		# cmd = subprocess.Popen(['docker','exec','-it','gatk','gatk','MarkDuplicatesSpark','--verbosity','ERROR','-I','%s/%s.sam' % (prep_folder,fastq_data[fastqfile]['sample']),'-O','%s/%s.markdup.sorted.bam' % (prep_folder,fastq_data[fastqfile]['sample'])], stdout=open('%s/markduplicates.stdout.txt' % prep_folder,'w'), stderr=open('%s/markduplicates.stderr.txt' % prep_folder,'w'))
		# cmd.communicate()

		# QUALITY SCORE RECALIBRATION (KEEP KNOWN SITES?)
		logging.info("\t\t- [%s] Base recalibration ..." % time.strftime("%H:%M:%S"))
		# cmd = subprocess.Popen(['docker','exec','-it','gatk','gatk','BaseRecalibratorSpark','-I','%s/%s.markdup.sorted.bam' % (prep_folder,fastq_data[fastqfile]['sample']),'-R',fastq_data[fastqfile]['reference'],'--known-sites','%s/reference_files/mutect/af-only-gnomad.raw.sites.b37.vcf.gz' % pipeline_folder,'-O','%s/%s.recal_data.table' % (prep_folder,fastq_data[fastqfile]['sample'])], stdout=open('%s/baserecalibrator.stdout.txt' % prep_folder,'w'), stderr=open('%s/baserecalibrator.stderr.txt' % prep_folder,'w'))
		# cmd.communicate()
		# cmd = subprocess.Popen(['docker','exec','-it','gatk','gatk','ApplyBQSRSpark','-I','%s/%s.markdup.sorted.bam' % (prep_folder,fastq_data[fastqfile]['sample']),'-R',fastq_data[fastqfile]['reference'],'--bqsr-recal-file','%s/%s.recal_data.table' % (prep_folder,fastq_data[fastqfile]['sample']),'-O',bam], stdout=open('%s/applybqsr.stdout.txt' % prep_folder,'w'), stderr=open('%s/applybqsr.stderr.txt' % prep_folder,'w'))
		# cmd.communicate()
	else:
		if paired_end_mbc:
			# PROCESS MBC
			logging.info("\t\t- [%s] Process MBC ..." % time.strftime("%H:%M:%S"))
			# cmd = subprocess.Popen(['%s/agent/agent.sh' % pipeline_folder,'locatit','-i','-C','-m','1','-U','-r','-IS','-OB','-c','100','-q','20','-c','100','-l',fastq_data[fastqfile]['covered_bed'],'-o','%s/%s.markdup.sorted.bam' % (prep_folder,fastq_data[fastqfile]['sample']),'%s/%s.sam' % (prep_folder,fastq_data[fastqfile]['sample']),fastq_data[fastqfile]['fastq_r2']], stdout=open('%s/locatit.stdout.txt' % prep_folder,'w'), stderr=open('%s/locatit.stderr.txt' % prep_folder,'w'))
			# cmd.communicate()
		else:
			# SAM TO BAM
			logging.info("\t\t- [%s] Converting SAM to BAM ..." % time.strftime("%H:%M:%S"))
			# cmd = subprocess.Popen(['samtools','view','-@','12','-bS','-o','%s/%s.markdup.sorted.bam' % (prep_folder,fastq_data[fastqfile]['sample']),'%s/%s.sam' % (prep_folder,fastq_data[fastqfile]['sample'])], stdout=open('%s/samtools_view.stdout.txt' % prep_folder,'w'), stderr=open('%s/samtools_view.stderr.txt' % prep_folder,'w'))
			# cmd.communicate()

		# SORT BAM
		logging.info("\t\t- [%s] Sorting BAM ..." % time.strftime("%H:%M:%S"))
		# cmd = subprocess.Popen(['samtools','sort','-@','12','-O','bam','-o',bam,'%s/%s.markdup.sorted.bam' % (prep_folder,fastq_data[fastqfile]['sample'])], stdout=open('%s/samtools_sort.stdout.txt' % prep_folder,'w'), stderr=open('%s/samtools_sort.stderr.txt' % prep_folder,'w'))
		# cmd.communicate()

		# INDEX
		logging.info("\t\t- [%s] Indexing BAM ..." % time.strftime("%H:%M:%S"))
		# cmd = subprocess.Popen(['samtools','index',bam], stdout=open('%s/samtools_index.stdout.txt' % prep_folder,'w'), stderr=open('%s/samtools_index.stderr.txt' % prep_folder,'w'))
		# cmd.communicate()

	# remove SAM and BAM temp files (if bam correctly produced) # A METTRE PLUTOT EN FIN D'ANALYSE
	# if os.path.exists(bam):
		# subprocess.call(['rm','%s/%s.sam' % (prep_folder,fastq_data[fastqfile]['sample'])],stdout=FNULL,stderr=FNULL)
		# subprocess.call(['rm','%s/%s.markdup.sorted.bam' % (prep_folder,fastq_data[fastqfile]['sample'])],stdout=FNULL,stderr=FNULL)
		# subprocess.call(['rm','%s/%s.markdup.sorted.bam.bai' % (prep_folder,fastq_data[fastqfile]['sample'])],stdout=FNULL,stderr=FNULL)
		# subprocess.call(['rm','%s/%s.markdup.sorted.bam.sbi' % (prep_folder,fastq_data[fastqfile]['sample'])],stdout=FNULL,stderr=FNULL)

	fastq_data[fastqfile]['bam'] = bam

	###  __   __        ___  __        __   ___                             __     __  ###
	### /  ` /  \ \  / |__  |__)  /\  / _` |__      /\  |\ |  /\  |    \ / /__` | /__` ###
	### \__, \__/  \/  |___ |  \ /~~\ \__> |___    /~~\ | \| /~~\ |___  |  .__/ | .__/ ###

	logging.info("\t\t- [%s] coverageAnalysis ..." % (time.strftime("%H:%M:%S")))
	cmd_list = []

	coverage_folder = '%s/coverage' % fastq_data[fastqfile]['intermediate_folder']
	if not os.path.isdir(coverage_folder):
		subprocess.call(['mkdir', coverage_folder])

	# cmd = subprocess.Popen([
		# 'bash', '%s/coverageAnalysis/run_coverage_analysis.sh' % pipeline_folder,
		# '-L', 'hg19', # attention ABL1?
		# '-dg', # '-g' seulement ? a pour amplicons. # -d = ignore duplicates (pour target coverage)
		# '-D',coverage_folder,
		# '-B',fastq_data[fastqfile]['target_bed'],
		# fastq_data[fastqfile]['reference'],
		# fastq_data[fastqfile]['bam']
		# ],
		# stdout=open('%s/run_coverage_analysis.stdout.txt' % coverage_folder,'w'), 
		# stderr=open('%s/run_coverage_analysis.stderr.txt' % coverage_folder,'w'))
	# cmd.communicate()

	#also compute samtools depth
	# cmd = subprocess.Popen(['samtools', 'depth','-b',fastq_data[fastqfile]['target_bed'],fastq_data[fastqfile]['bam']],stdout=open('%s/depth.txt' % coverage_folder,'w'), stderr=open('%s/samtools_depth.stderr.txt' % coverage_folder,'w'))
	# cmd.communicate()

# # QUALIMAP
	# logging.info("\t\t- [%s] qualimap ..." % (time.strftime("%H:%M:%S")))
	# qualimap_folder = '%s/qualimap' % fastq_data[fastqfile]['intermediate_folder']
	# qualimap_sub_folder = '%s/%s_%s' % (qualimap_folder,fastq_data[fastqfile]['sample'],fastq_data[fastqfile]['barcode'])
	# if not os.path.isdir(qualimap_folder):
		# subprocess.call(['mkdir', qualimap_folder])
	# if not os.path.isdir(qualimap_sub_folder):
		# subprocess.call(['mkdir', qualimap_sub_folder])

	# cmd = subprocess.Popen([
		# 'bash', '%s/qualimap/qualimap' % pipeline_folder,
		# 'bamqc',
		# '-bam', fastq_data[fastqfile]['bam'],
		# '--genome-gc-distr','hg19',
		# '--skip-duplicated', 
		# '--collect-overlap-pairs',
		# '--feature-file', fastq_data[fastqfile]['target_bed'],
		# '-outdir', qualimap_sub_folder,
		# '--java-mem-size=2G'
		# ],
		# stdout=open('%s/qualimap.stdout.txt' % qualimap_sub_folder,'w'), 
		# stderr=open('%s/qualimap.stderr.txt' % qualimap_sub_folder,'w'))
# cmd.communicate()

# SAMTOOLS STATS
	logging.info("\t\t- [%s] samtools stats ..." % (time.strftime("%H:%M:%S")))
	cmd = subprocess.Popen([
		'samtools', 'stats',
		'-d',
		'-t', fastq_data[fastqfile]['target_bed'],#'-@', '2',
		fastq_data[fastqfile]['bam']
		],
		stdout=open('%s/%s_%s.stats.txt' % (fastq_data[fastqfile]['intermediate_folder'],fastq_data[fastqfile]['sample'],fastq_data[fastqfile]['barcode']),'w'))

# SAMTOOLS FLAGSTAT
	# logging.info("\t\t- [%s] samtools flagstat ..." % (time.strftime("%H:%M:%S")))
	# cmd = subprocess.Popen([
		# 'samtools', 'flagstat',
		# '-@', '2',
		# fastq_data[fastqfile]['bam']
		# ],
		# stdout=open('%s/%s_%s.flagstat.txt' % (fastq_data[fastqfile]['intermediate_folder'],fastq_data[fastqfile]['sample'],fastq_data[fastqfile]['barcode']),'w'))

# MOSDEPTH
	# logging.info("\t\t- [%s] mosdepth ..." % (time.strftime("%H:%M:%S")))
	# mosdepth_folder = '%s/mosdepth' % fastq_data[fastqfile]['intermediate_folder']
	# if not os.path.isdir(mosdepth_folder):
		# subprocess.call(['mkdir', mosdepth_folder])
	# os.chdir(mosdepth_folder)
	# cmd = subprocess.Popen([
		# '%s/mosdepth/mosdepth' % pipeline_folder,
		# '-b', fastq_data[fastqfile]['target_bed'],
		# '%s_%s' % (fastq_data[fastqfile]['sample'],fastq_data[fastqfile]['barcode']),
		# fastq_data[fastqfile]['bam']
		# ],
		# stdout=open('%s/mosdepth.stdout.txt' % mosdepth_folder,'w'), 
		# stderr=open('%s/mosdepth.stderr.txt' % mosdepth_folder,'w'))
# cmd.communicate()

# QC : FASTQC BAM
	# logging.info("\t\t- [%s] FastQC (BAM)..." % (time.strftime("%H:%M:%S")))
	# fastqc_folder = '%s/fastqc' % fastq_data[fastqfile]['intermediate_folder']
	# if not os.path.isdir(fastqc_folder):
		# subprocess.call(['mkdir', fastqc_folder])
	# cmd = subprocess.Popen(['perl','%s/FastQC/fastqc' % pipeline_folder,'--outdir',fastqc_folder,fastq_data[fastqfile]['bam']],stdout=open('%s/fastqc_bam.stdout.txt' % fastqc_folder,'w'),stderr=open('%s/fastqc_bam.stderr.txt' % fastqc_folder,'w'))

#####  __        __  ___  __   __        ___  __        __   ___ ###
##### |__) |    /  \  |  /  ` /  \ \  / |__  |__)  /\  / _` |__  ###
##### |    |___ \__/  |  \__, \__/  \/  |___ |  \ /~~\ \__> |___ ###

# if options.run and os.path.isfile('%s/barcodes.json' % run_folder):
	# logging.info(" [%s] plotCoverage ..." % time.strftime("%H:%M:%S"))
	# cmd = subprocess.Popen(['python', '%s/plotCoverage/plotCoverage.py' % pipeline_folder,'--run-folder', run_folder],stdout=open('%s/_plotCoverage/plotCoverage.stdout.txt' % run_folder,'w'), stderr=open('%s/_plotCoverage/plotCoverage.stderr.txt' % run_folder,'w'))
	# cmd.communicate()
	# with open('%s/_plotCoverage/plotCoverage.stderr.txt' % run_folder,'r') as stderr:
		# for line in stderr.readlines():
			# print line.replace('\n','')

#####  __            ###
##### /  ` |\ | \  / ###
##### \__, | \|  \/  ###

# if options.run:
	# logging.info(" [%s] CNV Analysis ..." % (time.strftime("%H:%M:%S")))
	# if not os.path.isdir('%s/_CNA' % run_folder):
		# subprocess.call(['mkdir', '%s/_CNA' % run_folder])
	# cmd = subprocess.Popen(['python','%s/CNV/run_cna.py' % pipeline_folder, run_folder], stdout=open('%s/_CNA/run_cna.stdout.txt' % run_folder,'w'), stderr=open('%s/_CNA/run_cna.stderr.txt' % run_folder,'w'))
	# cmd.communicate()
	# with open('%s/_CNA/run_cna.stderr.txt' % run_folder,'r') as stderr:
		# for line in stderr.readlines():
			# print line.replace('\n','')

# QUALIMAP
# logging.info("\n- [%s] Qualimap ..." % (time.strftime("%H:%M:%S")))
# qualimap_folder = '%s/_qualimap' % run_folder
# if not os.path.isdir(qualimap_folder):
	# subprocess.call(['mkdir', qualimap_folder])

# qualimap_data_path = '%s/qualimap-data.tsv' % qualimap_folder
# qualimap_input = open(qualimap_data_path,'w')
# for fastqfile in fastq_data:
	# qualimap_input.write('%s_%s\t%s\n' % (fastq_data[fastqfile]['sample'],fastq_data[fastqfile]['barcode'],fastq_data[fastqfile]['bam']))
# qualimap_input.close()

# cmd = subprocess.Popen([
	# 'bash', '%s/qualimap/qualimap' % pipeline_folder,
	# 'multi-bamqc',
	# '-data', qualimap_data_path,
	# '--genome-gc-distr','hg19',
	# '--skip-duplicated', #--collect-overlap-pairs??
	# '--feature-file', fastq_data[fastqfile]['target_bed'],
	# '-outdir', qualimap_folder,
	# '--java-mem-size=12G'
	# ],
	# stdout=open('%s/qualimap.stdout.txt' % qualimap_folder,'w'), 
	# stderr=open('%s/qualimap.stderr.txt' % qualimap_folder,'w'))
# cmd.communicate()

##### #####            __              ___     __                         __                __                     __  ___      ___    __       ##### ##### 
##### ##### \  /  /\  |__) |  /\  |\ |  |  __ /  `  /\  |    |    | |\ | / _`     /\  |\ | |  \     /\  |\ | |\ | /  \  |   /\   |  | /  \ |\ | ##### ##### 
##### #####  \/  /~~\ |  \ | /~~\ | \|  |     \__, /~~\ |___ |___ | | \| \__>    /~~\ | \| |__/    /~~\ | \| | \| \__/  |  /~~\  |  | \__/ | \| ##### ##### 
##### #####                                                                                                                                     ##### ##### 

exit()
logging.info("\n- [%s] VARIANT-CALLING AND ANNOTATION ..." % (time.strftime("%H:%M:%S")))
for fastqfile in fastq_data:
	logging.info("\t- %s :" % (fastq_data[fastqfile]['sample']))

	# FOR SKIPING SAMPLES # WARNING check-contamination will not work #
	if sampleask:
		proceed = raw_input('continue? (y/n/stopask)\n')              
		if proceed == 'y' or proceed == 'yes':                        
			pass
		elif proceed == 'stopask':
			sampleask = False
			pass                                                       
		else:                                                         
			continue                                                  

	### __   ___  ___  __             __              ___ 
	### |  \ |__  |__  |__) \  /  /\  |__) |  /\  |\ |  |  
	### |__/ |___ |___ |     \/  /~~\ |  \ | /~~\ | \|  |  

	logging.info("\t\t - [%s] DeepVariant ..." % time.strftime("%H:%M:%S"))
	deepvariant_folder = '%s/deepvariant' % fastq_data[fastqfile]['intermediate_folder']
	if not os.path.isdir(deepvariant_folder):
		subprocess.call(['mkdir', deepvariant_folder])

	cmd = subprocess.Popen(['docker','exec','-it','deepvariant','/opt/deepvariant/bin/run_deepvariant',
	'--model_type=WES',
	'--ref=%s' % fastq_data[fastqfile]['reference'],
	'--reads=%s' % fastq_data[fastqfile]['bam'],
	'--output_vcf=%s/%s.deepvariant.vcf' % (deepvariant_folder,fastq_data[fastqfile]['sample']),
	'--regions=%s' % fastq_data[fastqfile]['merged_bed'],
	'--num_shards=12'
	],stdout=open('%s/deepvariant.stdout.txt'%deepvariant_folder,'w'),stderr=open('%s/deepvariant.stderr.txt'%deepvariant_folder,'w'))
	cmd.communicate()

	###            __   __   __            ###
	### \  /  /\  |__) /__` /  `  /\  |\ | ###
	###  \/  /~~\ |  \ .__/ \__, /~~\ | \| ###

	logging.info("\t\t - [%s] VarScan2 ..." % time.strftime("%H:%M:%S"))
	varscan2_folder = '%s/varscan2' % fastq_data[fastqfile]['intermediate_folder']
	if not os.path.isdir(varscan2_folder):
		subprocess.call(['mkdir', varscan2_folder])

	logging.info("\t\t\t - [%s] mpileup ..." % time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen(['samtools','mpileup','-f',fastq_data[fastqfile]['reference'],'-l',fastq_data[fastqfile]['target_bed'],'-o','%s/%s.pileup' % (varscan2_folder,fastq_data[fastqfile]['sample']),fastq_data[fastqfile]['bam']],stdout=open('%s/mpileup.stdout.txt'%varscan2_folder,'w'),stderr=open('%s/mpileup.stderr.txt'%varscan2_folder,'w'))
	cmd.communicate()

	logging.info("\t\t\t - [%s] mpileup2cns ..." % time.strftime("%H:%M:%S"))
	cmd1 = subprocess.Popen(['java','-jar','%s/varscan2/VarScan.v2.4.2.jar' % pipeline_folder,'mpileup2cns','%s/%s.pileup' % (varscan2_folder,fastq_data[fastqfile]['sample']),'--min-var-freq','0.02','--output-vcf','1','--variants'],stdout=open('%s/%s.varscan2.cns.vcf' % (varscan2_folder,fastq_data[fastqfile]['sample']),'w'),stderr=open('%s/varscan2.cns.stderr.txt'%varscan2_folder,'w'))
	cmd1.communicate()

	logging.info("\t\t\t - [%s] filter ..." % time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen(['java','-jar','%s/varscan2/VarScan.v2.4.2.jar' % pipeline_folder,'filter','%s/%s.varscan2.cns.vcf' % (varscan2_folder,fastq_data[fastqfile]['sample']),
	'--min-coverage','10',
	'--min-reads2', '2',
	'--min-strands2','1',
	'--min-avg-qual','10',
	'--min-var-freq','0.01',
	'--output-file','%s/%s.varscan2.filtered.vcf' % (varscan2_folder,fastq_data[fastqfile]['sample'])
	],stdout=open('%s/varscan2.filter.stdout.txt'%varscan2_folder,'w'),stderr=open('%s/varscan2.filter.stderr.txt'%varscan2_folder,'w'))
	cmd.communicate()

	###       __   ___  __   ___  __  ###
	### |    /  \ |__  |__) |__  /  \ ###
	### |___ \__/ |    |  \ |___ \__X ###

	lofreq_folder = '%s/lofreq' % fastq_data[fastqfile]['intermediate_folder']
	if not os.path.isdir(lofreq_folder):
		subprocess.call(['mkdir', lofreq_folder])

	logging.info("\t\t - [%s] LoFreq ..." % time.strftime("%H:%M:%S"))

	# LOFREQ INDELQUAL dindel
	logging.info("\t\t\t - [%s] indelqual ..." % time.strftime("%H:%M:%S"))
	if fastq_data[fastqfile]['fastq_r2'] != None:
		cmd = subprocess.Popen(['%s/lofreq/bin/lofreq' % pipeline_folder,'indelqual','--dindel','-f',fastq_data[fastqfile]['reference'],'-o','%s/%s.indelqual.bam' % (lofreq_folder,fastq_data[fastqfile]['sample']),fastq_data[fastqfile]['bam']], stdout=open('%s/lofreq_indelqual_dindel.stdout.txt' % lofreq_folder,'w'), stderr=open('%s/lofreq_indelqual_dindel.stderr.txt' % lofreq_folder,'w'))
		cmd.communicate()
	else: # Set INDELQUAL uniform 45 for ion torrent
		cmd = subprocess.Popen(['%s/lofreq/bin/lofreq' % pipeline_folder,'indelqual','--uniform','45','-o','%s/%s.indelqual.bam' % (lofreq_folder,fastq_data[fastqfile]['sample']),fastq_data[fastqfile]['bam']], stdout=open('%s/lofreq_indelqual_uniform.stdout.txt' % lofreq_folder,'w'), stderr=open('%s/lofreq_indelqual_uniform.stderr.txt' % lofreq_folder,'w'))
		cmd.communicate()
	cmd = subprocess.Popen(['samtools','index','%s/%s.indelqual.bam' % (lofreq_folder,fastq_data[fastqfile]['sample'])], stdout=open('%s/samtools_index.stdout.txt' % lofreq_folder,'w'), stderr=open('%s/samtools_index.stderr.txt' % lofreq_folder,'w'))
	cmd.communicate()
 
	logging.info("\t\t\t - [%s] call ..." % time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen(['%s/lofreq/bin/lofreq' % pipeline_folder,'call-parallel','--pp-threads','12','--call-indels','-f',fastq_data[fastqfile]['reference'],'-l',fastq_data[fastqfile]['target_bed'],'-o','%s/%s.lofreq.vcf' % (lofreq_folder,fastq_data[fastqfile]['sample']),'%s/%s.indelqual.bam' % (lofreq_folder,fastq_data[fastqfile]['sample'])],stdout=open('%s/lofreq.stdout.txt'%lofreq_folder,'w'),stderr=open('%s/lofreq.stderr.txt'%lofreq_folder,'w'))
	cmd.communicate()

	logging.info("\t\t\t - [%s] filter ..." % time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen(['%s/lofreq/bin/lofreq' % pipeline_folder,'filter',
	'-i','%s/%s.lofreq.vcf' % (lofreq_folder,fastq_data[fastqfile]['sample']),
	'-o','%s/%s.lofreq.filtered.vcf' % (lofreq_folder,fastq_data[fastqfile]['sample']),
	'--cov-min','10',
	'--af-min','0.01',
	'--snvqual-thresh','10'
	],stdout=open('%s/lofreq.filter.stdout.txt'%lofreq_folder,'w'),stderr=open('%s/lofreq.filter.stderr.txt'%lofreq_folder,'w'))
	cmd.communicate()

	# remove temp indelqual bam
	subprocess.call(['rm','%s/%s.indelqual.bam' % (lofreq_folder,fastq_data[fastqfile]['sample'])])

	###           ___  ___  __  ___ ###
	### |\/| |  |  |  |__  /  `  |  ###
	### |  | \__/  |  |___ \__,  |  ###

	logging.info("\t\t - [%s] Mutect2 ..." % time.strftime("%H:%M:%S"))
	mutect2_folder = '%s/mutect2' % fastq_data[fastqfile]['intermediate_folder']
	if not os.path.isdir(mutect2_folder):
		subprocess.call(['mkdir', mutect2_folder])

	do_chunk = True
	mutect_chunk = 10
	if do_chunk:
		logging.info("\t\t\t - [%s] SplitIntervals ..." % time.strftime("%H:%M:%S"))
		cmd = subprocess.Popen(['docker','exec','-it','gatk','gatk','SplitIntervals','-R',fastq_data[fastqfile]['reference'],'-L',fastq_data[fastqfile]['intervals'],'-O',mutect2_folder,'--scatter-count',str(mutect_chunk)],stdout=open('%s/splitintervals.stdout.txt'%mutect2_folder,'w'),stderr=open('%s/splitintervals.stderr.txt'%mutect2_folder,'w'))
		cmd.communicate()

		ps = []
		vcf_chunk_list = []
		logging.info("\t\t\t - [%s] Running %s chunks ..." % (time.strftime("%H:%M:%S"),mutect_chunk))
		FNULL = open(os.devnull, 'w')
		for i in range(mutect_chunk):
			# interval_chunk = fastq_data[fastqfile]['intervals'].replace('.intervals','.%04d-scattered.intervals'%i)
			interval_chunk = '%s/%04d-scattered.interval_list'%(mutect2_folder,i)
			vcf_chunk = '%s/%s.mutect2.chunk%04d.vcf' % (mutect2_folder,fastq_data[fastqfile]['sample'],i)
			vcf_chunk_list.append(vcf_chunk)
			cmd = subprocess.Popen(['docker','exec','-it','gatk','gatk','Mutect2',
			'-I',fastq_data[fastqfile]['bam'],
			'-R',fastq_data[fastqfile]['reference'],
			'-L',interval_chunk,
			'--max-mnp-distance','0',
			'--max-reads-per-alignment-start','0',
			'-O',vcf_chunk],stdout=open('%s/mutect2.chunk%04d.stdout.txt' % (mutect2_folder,i),'w'),stderr=open('%s/mutect2.chunk%04d.stderr.txt' % (mutect2_folder,i),'w'))
			ps.append(cmd)
		for p in ps:
			p.wait()
		subprocess.call(['stty','sane'])

		logging.info("\t\t\t - [%s] GatherVcfs ..." % time.strftime("%H:%M:%S"))
		args = ['docker','exec','-it','gatk','gatk','GatherVcfs','-O','%s/%s.mutect2.vcf' % (mutect2_folder,fastq_data[fastqfile]['sample'])]
		for vcf_chunk in vcf_chunk_list:
			args.append('-I')
			args.append(vcf_chunk)
		cmd = subprocess.Popen(args,stdout=open('%s/gathervcf.stdout.txt' % mutect2_folder,'w'),stderr=open('%s/gathervcf.stderr.txt' % mutect2_folder,'w'))
		cmd.wait()

		logging.info("\t\t\t - [%s] MergeMutectStats ..." % time.strftime("%H:%M:%S"))
		args = ['docker','exec','-it','gatk','gatk','MergeMutectStats','-O','%s/%s.mutect2.vcf.stats' % (mutect2_folder,fastq_data[fastqfile]['sample'])]
		for vcf_chunk in vcf_chunk_list:
			args.append('-stats')
			args.append(vcf_chunk+'.stats')
		cmd = subprocess.Popen(args,stdout=open('%s/mergemutectstats.stdout.txt' % mutect2_folder,'w'),stderr=open('%s/mergemutectstats.stderr.txt' % mutect2_folder,'w'))
		cmd.wait()
	else:
		cmd = subprocess.Popen('docker exec -it gatk gatk Mutect2 -I %s -R %s -L %s --max-mnp-distance 0 --max-reads-per-alignment-start 0 --min-base-quality-score 1 -O %s/%s.mutect2.vcf' % (fastq_data[fastqfile]['bam'],fastq_data[fastqfile]['reference'],fastq_data[fastqfile]['intervals'],mutect2_folder,fastq_data[fastqfile]['sample']),shell=True,stdout=open('%s/mutect2.stdout.txt'%mutect2_folder,'w'),stderr=open('%s/mutect2.stderr.txt'%mutect2_folder,'w'))
		cmd.communicate()

	#print "\t\t - [%s] Tabix ..." % time.strftime("%H:%M:%S")
	#subprocess.call(['tabix','-p','vcf',mutect2],stdout=open('%s/tabix.stdout.txt' % mutect2_folder,'w'),stderr=open('%s/tabix.stderr.txt' % mutect2_folder,'w'))

	logging.info("\t\t\t - [%s] FilterMutectCalls ..." % time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen([
		'docker','exec','-it','gatk','gatk', 'FilterMutectCalls',
		'-V','%s/%s.mutect2.vcf' % (mutect2_folder,fastq_data[fastqfile]['sample']),
		'-R',fastq_data[fastqfile]['reference'],
		'--max-events-in-region', '20',
		'--min-allele-fraction', '0.01',
		'--min-median-base-quality','10',
		'--unique-alt-read-count','10',
		'-O','%s/%s.mutect2.filtered.vcf' % (mutect2_folder,fastq_data[fastqfile]['sample'])
		#'--panel-of-normals', '%s/reference_files/mutect/colon_lung_pon.vcf.gz' % pipeline_folder
	], stdout=open('%s/filtermutectcalls.stdout.txt' % mutect2_folder,'w'),stderr=open('%s/filtermutectcalls.stderr.txt' % mutect2_folder,'w'))
	cmd.communicate()

###         __   ___  __  ___     __   __  ### 
### | |\ | /__` |__  |__)  |     |  \ |__) ### 
### | | \| .__/ |___ |  \  |     |__/ |__) ### 

	# add variantmetrics for each variant (if analysis already done, delete old variantmetrics first)
	# add new entries for new variants in db
	logging.info("\t\t - [%s] insert new variants and metrics into DB ..." % (time.strftime("%H:%M:%S")))
	abl1 = 'no'
	if 'ABL1_NM_005157.fasta' in fastq_data[fastqfile]['reference']:
		abl1 = 'yes'

	cmd = subprocess.Popen(['python','%s/variantBase/insert_db_variants_illumina.py' % pipeline_folder,
		'--analysis', fastq_data[fastqfile]['analysis_id'],
		'--vcf', '%s/%s.mutect2.filtered.vcf' % (mutect2_folder,fastq_data[fastqfile]['sample']),
		'--vcf', '%s/%s.varscan2.filtered.vcf' % (varscan2_folder,fastq_data[fastqfile]['sample']),
		'--vcf', '%s/%s.lofreq.filtered.vcf' % (lofreq_folder,fastq_data[fastqfile]['sample']),
		'--vcf', '%s/%s.deepvariant.vcf' % (deepvariant_folder,fastq_data[fastqfile]['sample']),
		'--abl1', abl1
		], 
		stdout=open('%s/insert_db_variants.stdout.txt' % fastq_data[fastqfile]['intermediate_folder'],'w'), 
		stderr=open('%s/insert_db_variants.stderr.txt' % fastq_data[fastqfile]['intermediate_folder'],'w'))
	cmd.communicate()
	with open('%s/insert_db_variants.stdout.txt' % fastq_data[fastqfile]['intermediate_folder'],'r') as stdout:
		for line in stdout.readlines():
			print '\t\t\t' + line.replace('\n','')
	with open('%s/insert_db_variants.stderr.txt' % fastq_data[fastqfile]['intermediate_folder'],'r') as stderr:
		for line in stderr.readlines():
			print '\t\t\t' + line.replace('\n','')

	###                 __  ___      ___  ___          ___                    __              ___  __  ###
	###  /\  |\ | |\ | /  \  |   /\   |  |__     |\ | |__  |  |    \  /  /\  |__) |  /\  |\ |  |  /__` ###
	### /~~\ | \| | \| \__/  |  /~~\  |  |___    | \| |___ |/\|     \/  /~~\ |  \ | /~~\ | \|  |  .__/ ###

	logging.info("\t\t - [%s] annotate new DB variants ..." % (time.strftime("%H:%M:%S")))

	logging.info("\t\t\t - [%s] HGVS check ..." % (time.strftime("%H:%M:%S")))
	cmd = subprocess.Popen(['python','%s/variantAnnotation/annotate_variantbase.step0.py' % pipeline_folder,'--new'], #--new only annote new variants
	stdout=open('%s/annotate_variantbase.step0.stdout.txt' % fastq_data[fastqfile]['intermediate_folder'],'w'), 
	stderr=open('%s/annotate_variantbase.step0.stderr.txt' % fastq_data[fastqfile]['intermediate_folder'],'w'))
	cmd.communicate()
	with open('%s/annotate_variantbase.step0.stderr.txt' % fastq_data[fastqfile]['intermediate_folder'],'r') as stderr:
		for line in stderr.readlines():
			if not line.startswith('No handlers could be found for logger "hgvs"'):
				print '\t\t\t' + line.replace('\n','')

	logging.info("\t\t\t - [%s] Annovar and VEP ..." % (time.strftime("%H:%M:%S")))
	cmd = subprocess.Popen(['python','%s/variantAnnotation/annotate_variantbase.step1.py' % pipeline_folder,'--new','--output-folder',fastq_data[fastqfile]['intermediate_folder']], 
	stdout=open('%s/annotate_variantbase.step1.stdout.txt' % fastq_data[fastqfile]['intermediate_folder'],'w'), 
	stderr=open('%s/annotate_variantbase.step1.stderr.txt' % fastq_data[fastqfile]['intermediate_folder'],'w'))
	cmd.communicate()
	with open('%s/annotate_variantbase.step1.stderr.txt' % fastq_data[fastqfile]['intermediate_folder'],'r') as stderr:
		for line in stderr.readlines():
			print '\t\t\t' + line.replace('\n','')

	logging.info("\t\t\t - [%s] Merging annotations and updating DB ..." % (time.strftime("%H:%M:%S")))
	if os.path.isfile('%s/annovar/annovar_input.tsv' % fastq_data[fastqfile]['intermediate_folder']):
		cmd = subprocess.Popen(['python','%s/variantAnnotation/annotate_variantbase.step2.py' % pipeline_folder,'--annovar-results', '%s/annovar/annovar.hg19_multianno.txt' % fastq_data[fastqfile]['intermediate_folder'],'--vep-results', '%s/vep/vep_output.tsv' % fastq_data[fastqfile]['intermediate_folder']], 
		stdout=open('%s/annotate_variantbase.step2.stdout.txt' % fastq_data[fastqfile]['intermediate_folder'],'w'),
		stderr=open('%s/annotate_variantbase.step2.stderr.txt' % fastq_data[fastqfile]['intermediate_folder'],'w'))
		cmd.communicate()
		with open('%s/annotate_variantbase.step2.stderr.txt' % fastq_data[fastqfile]['intermediate_folder'],'r') as stderr:
			for line in stderr.readlines():
				print '\t\t\t' + line.replace('\n','')
				
###  ___                   __   ___  __   __   __  ___ ###
### |__  | |\ |  /\  |    |__) |__  |__) /  \ |__)  |  ###
### |    | | \| /~~\ |___ |  \ |___ |    \__/ |  \  |  ###

	if not os.path.isdir('%s/finalReport' % fastq_data[fastqfile]['intermediate_folder']):
		subprocess.call(['mkdir', '%s/finalReport' % fastq_data[fastqfile]['intermediate_folder']])
	finalReport_path = '%s/%s_%s.finalReport.xlsx' % (fastq_data[fastqfile]['sample_folder'],fastq_data[fastqfile]['sample'],fastq_data[fastqfile]['barcode'])
	
	xmin = '300' # for coverage analysis sheet in finalreport
	if 'cDNA' in param_hotspot_only:
		xmin = '2000'
		
	logging.info("\t\t - [%s] finalReport ..." % (time.strftime("%H:%M:%S")))
	cmd = subprocess.Popen([
		'python','%s/finalReport/finalReport.py' % pipeline_folder,
		'--analysis', fastq_data[fastqfile]['analysis_id'],
		'--xmin', xmin,
	], 
	stdout=open('%s/finalReport/finalReport.stdout.txt' % fastq_data[fastqfile]['intermediate_folder'],'w'), 
	stderr=open('%s/finalReport/finalReport.stderr.txt' % fastq_data[fastqfile]['intermediate_folder'],'w'))
	cmd.communicate()
	with open('%s/finalReport/finalReport.stdout.txt' % fastq_data[fastqfile]['intermediate_folder'],'r') as stdout:
		for line in stdout.readlines():
			print '\t\t\t' + line.replace('\n','')
	with open('%s/finalReport/finalReport.stderr.txt' % fastq_data[fastqfile]['intermediate_folder'],'r') as stderr:
		for line in stderr.readlines():
			print '\t\t\t' + line.replace('\n','')

	#################################################

# MUTLIQC
if options.run:
	logging.info("\t- [%s] MultiQC ..." % (time.strftime("%H:%M:%S")))
	multiqc_folder = '%s/_multiqc' % run_folder
	if not os.path.isdir(multiqc_folder):
		subprocess.call(['mkdir', multiqc_folder])
	cmd = subprocess.Popen(['multiqc',run_folder,'--outdir',multiqc_folder],stdout=open('%s/multiqc.stdout.txt' % multiqc_folder,'w'),stderr=open('%s/multiqc.stderr.txt' % multiqc_folder,'w'))
	cmd.communicate()

