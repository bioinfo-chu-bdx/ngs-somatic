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


def check_stderr(stderr_path,indent=0):
	with open(stderr_path,'r') as stderr:
		for line in stderr.readlines():
			line = line.replace('\n','')
			for i in range(indent):
				line = '\t%s' % line
			logging.info(line)


### GATHERING PARAMETERS ############################################################

parser = OptionParser()
parser.add_option('-b', '--bam', help="Input bam file for SINGLE BAM ANALYSIS (must contain flow signal from TSS)", dest='bam')
parser.add_option('-r', '--run', help="Run folder path  for FULL RUN ANALYSIS", dest='run')
(options, args) = parser.parse_args()

if options.bam and options.run:
	sys.stderr.write("[run_analysis.py] Error: <--bam> and <--full-run> are not compatibles\n")
	sys.exit()
if options.run:
	bamlist = glob.glob('%s/*/*.bam' % options.run)
	bamlist = [item for item in bamlist if not 'processed' in item]
	run_folder = options.run
elif options.bam:
	bamlist = [options.bam]
	run_folder = os.path.dirname(os.path.dirname(options.bam))
else:
	sys.stderr.write("[run_analysis.py] Error: no <--bam> or <--full-run> specified\n")
	sys.exit()
	
# set up logging to file
if options.run:
	logging.basicConfig(level=logging.INFO,format='%(levelname)-7s %(message)s',filename='/%s/run_analysis.log' % run_folder,filemode='w')
else:
	logging.basicConfig(level=logging.INFO,format='%(levelname)-7s %(message)s',filename='/%s/bam_analysis.log' % os.path.dirname(options.bam),filemode='w')
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
bam_data = {}

if run_folder.endswith('/'):
	run_name = os.path.basename(os.path.dirname(run_folder))
else:
	run_name = os.path.basename(run_folder)
	
######################################################################################

logging.info(" [%s] Starting analysis ..." % (time.strftime("%H:%M:%S")))

# CREATE RUN DB ENTRY - IF INEXISTANT -
if options.run:
	db_cur.execute("SELECT runID FROM Run WHERE runID='%s'"%run_name)
	if db_cur.fetchone() is None:
		logging.info("\t - [%s] Create new Run DB entry ..." % (time.strftime("%H:%M:%S")))
		platform = 'S5' # default, a ameliorer
		if 'S5' in run_name:
			platform = 'S5'
		if 'PGM' in run_name:
			platform = 'PGM'
		if os.path.isfile('%s/pgm_logs.zip' % run_folder):
			archive = zipfile.ZipFile('%s/pgm_logs.zip' % run_folder, 'r')
			with archive.open('explog_final.txt') as explog_final:
				for line in explog_final:
					if line.startswith('Start Time:'):
						start_time = line.split('Start Time:')[-1].strip()
						if platform == 'PGM': # format = Tue Jan 24 08:54:47 2017
							start_time = start_time.replace('  ',' ')
							if start_time.startswith(' '):
								start_time = start_time[1:]
							start_time = start_time.split(' ')
							run_date = date(int(start_time[4]),month2num[start_time[1]],int(start_time[2]))
						elif platform == 'S5': # format = 01/18/2017 10:23:51
							start_time = start_time.replace(' ','/').split('/')
							run_date = date(int(start_time[2]),int(start_time[0]),int(start_time[1]))
						break
		else:
			run_date = date.today()
		try:
			db_cur.execute("INSERT INTO Run (runID, platform, runPath, runDate) VALUES ('%s', '%s', '%s', '%s')" % (run_name,platform,run_folder,run_date))
			db_con.commit()
		except Exception as e:
			logging.warning("*WARNING* (RUN table)** %s" % e)
			sys.exit()
		
### RUN ANALYSIS ######################################################################		

for bamfile in sorted(bamlist) :
	sample = bamfile.split('/')[-1].split('_IonXpress')[0]
	barcode = 'IonXpress_%s' % bamfile.split('IonXpress_')[-1].split('.bam')[0]
	sample_folder = os.path.dirname(bamfile)
	
	## CREATE SAMPLE DB ENTRY - IF INEXISTANT -
	sample_id = barcodes_json[barcode]['sample_id']
	name = sample.replace(sample_id,'')
	while name.endswith('-') or name.endswith('_'):
		name = name[:-1]
	iscontrol = 0
	if sample_id.startswith('CONTROL-'):
		iscontrol = 1
	db_cur.execute("SELECT sampleID FROM Sample WHERE sampleID='%s'" % sample_id)
	if db_cur.fetchone() is None:
		logging.info("\t - [%s] new Sample and Analysis DB entry for %s ..." % (time.strftime("%H:%M:%S"), sample))
		try:
			db_cur.execute("INSERT INTO Sample (sampleID, sampleName, isControl) VALUES ('%s', '%s', '%s')" % (sample_id,name,iscontrol))
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
			db_cur.execute("INSERT INTO Analysis (analysisID, sample, barcode, run, panel, bamPath, analysisDate) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s')" % (analysis_id, sample_id, barcode, run_name, bed_name, bamfile, time.strftime("%Y-%m-%d")))
			db_con.commit()
		except Exception as e:
			logging.warning("\t*WARNING* (Analysis table)** %s" % e)
	else:	
		analysis_id = db_analysis['analysisID']
	logging.info("\t - %s, analysisID = %s" % (sample,analysis_id))
		
	### RUN TYPE PARAMETERS ###	
	if not run_type:
		logging.warning("\t -- Warning : run type (panel) not found for %s. Sample will not be processed." % sample)
		continue
	reference = global_param['run_type'][run_type]['reference']
	target_bed = global_param['run_type'][run_type]['target_bed']
	bed_merged = global_param['run_type'][run_type]['merged_bed']
	param = global_param['run_type'][run_type]['vc_parameters']
	param_hotspot_only = global_param['run_type'][run_type]['vc_parameters_hotspot_only']
	if 'PL' in sample_id : # cDNA sample
		if 'vc_parameters_hotspot_cdna' in global_param['run_type'][run_type]:
			logging.info("- Sample %s is cDNA" % sample)
			param_hotspot_only = global_param['run_type'][run_type]['vc_parameters_hotspot_cdna']
	hotspot_vcf = global_param['run_type'][run_type]['hotspot_vcf']
	checkconta_read_len = global_param['run_type'][run_type]['checkContamination_read_len']

	# DOSSIER INTERMEDIAIRE
	intermediate_folder = '%s/intermediate_files' % sample_folder
	if os.path.isfile('%s.zip' % intermediate_folder): # UNZIP ALL intermediates_files.zip (if RERUN)
		fnull = open(os.devnull, 'w')
		subprocess.call(['unzip', '-o', '%s.zip' % intermediate_folder,'-d', intermediate_folder],stdout=fnull)
	elif not os.path.isdir(intermediate_folder):
		subprocess.call(['mkdir', intermediate_folder])
		
	bam_data[bamfile] = {'barcode':barcode,'sample':sample,'sample_id':sample_id,'analysis_id':analysis_id,'run_type':run_type,'reference':reference,'sample_folder':sample_folder,'intermediate_folder':intermediate_folder,'target_bed':target_bed,'bed_merged':bed_merged,'param':param,'param_hotspot_only':param_hotspot_only,'hotspot_vcf':hotspot_vcf,'checkconta_read_len':checkconta_read_len}

db_con.close()
		
###  __   __        ___  __        __   ___                             __     __  ###
### /  ` /  \ \  / |__  |__)  /\  / _` |__      /\  |\ |  /\  |    \ / /__` | /__` ###
### \__, \__/  \/  |___ |  \ /~~\ \__> |___    /~~\ | \| /~~\ |___  |  .__/ | .__/ ###

#logging.info(" [%s] coverageAnalysis ..." % (time.strftime("%H:%M:%S")))
#cmd_list = []

#for bamfile in bam_data:
	#coverage_folder = '%s/coverage' % bam_data[bamfile]['intermediate_folder']
	#if not os.path.isdir(coverage_folder):
		#subprocess.call(['mkdir', coverage_folder])
	#cmd = subprocess.Popen([
		#'bash', '%s/coverageAnalysis/run_coverage_analysis.sh' % pipeline_folder,
		#'-L', 'hg19', # attention ABL1?
		#'-ag',
		#'-D',coverage_folder,
		#'-B',bam_data[bamfile]['target_bed'],
		#bam_data[bamfile]['reference'],
		#bamfile
		#], 	
		#stdout=open('%s/run_coverage_analysis.stdout.txt' % coverage_folder,'w'), 
		#stderr=open('%s/run_coverage_analysis.stderr.txt' % coverage_folder,'w'))
	#cmd_list.append((cmd,'%s/run_coverage_analysis.stderr.txt' % coverage_folder))
#for process in cmd_list:
	#process[0].communicate()
	#stderr_file = process[1]
	#check_stderr(stderr_file,indent=1)
	
######  __        __  ___  __   __        ___  __        __   ___ ###
###### |__) |    /  \  |  /  ` /  \ \  / |__  |__)  /\  / _` |__  ###
###### |    |___ \__/  |  \__, \__/  \/  |___ |  \ /~~\ \__> |___ ###

#if options.run and os.path.isfile('%s/barcodes.json' % run_folder):
	#logging.info(" [%s] plotCoverage ..." % time.strftime("%H:%M:%S"))
	#if not os.path.isdir('%s/_plotCoverage' % run_folder):
		#subprocess.call(['mkdir', '%s/_plotCoverage' % run_folder])
	#cmd = subprocess.Popen(['python', '%s/plotCoverage/plotCoverage.py' % pipeline_folder,'--run-folder', run_folder],stdout=open('%s/_plotCoverage/plotCoverage.stdout.txt' % run_folder,'w'), stderr=open('%s/_plotCoverage/plotCoverage.stderr.txt' % run_folder,'w'))
	#cmd.communicate()
	#check_stderr('%s/_plotCoverage/plotCoverage.stderr.txt'%run_folder,indent=1)

######  __            ###
###### /  ` |\ | \  / ###
###### \__, | \|  \/  ###

#if options.run:
	#logging.info(" [%s] CNV Analysis ..." % (time.strftime("%H:%M:%S")))
	#if not os.path.isdir('%s/_CNA' % run_folder):
		#subprocess.call(['mkdir', '%s/_CNA' % run_folder])
	#cmd = subprocess.Popen(['python','%s/CNV/run_cna.py' % pipeline_folder, run_folder], stdout=open('%s/_CNA/run_cna.stdout.txt' % run_folder,'w'), stderr=open('%s/_CNA/run_cna.stderr.txt' % run_folder,'w'))
	#cmd.communicate()
	#check_stderr('%s/_CNA/run_cna.stderr.txt'%run_folder,indent=1)

### RESTE DE L'ANALYSE PATIENT PAR PATIENT ###	
for bamfile in sorted(bamlist) :
	
	###################################################################
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
	###################################################################
	
	barcode = bam_data[bamfile]['barcode']
	sample = bam_data[bamfile]['sample']
	sample_id = bam_data[bamfile]['sample_id']
	analysis_id = bam_data[bamfile]['analysis_id']
	run_type = bam_data[bamfile]['run_type']
	reference = bam_data[bamfile]['reference']
	sample_folder = bam_data[bamfile]['sample_folder']
	intermediate_folder = bam_data[bamfile]['intermediate_folder']
	target_bed = bam_data[bamfile]['target_bed']
	bed_merged = bam_data[bamfile]['bed_merged']
	param = bam_data[bamfile]['param']
	param_hotspot_only = bam_data[bamfile]['param_hotspot_only']
	hotspot_vcf = bam_data[bamfile]['hotspot_vcf']
	checkconta_read_len = bam_data[bamfile]['checkconta_read_len']
	
	logging.info(" [%s] Processing %s - %s :" % (time.strftime("%H:%M:%S"),barcode,sample))
	logging.info("\t * panel %s" % run_type)
	logging.info("\t * analysisID %s" % analysis_id)
	
###            __              ___     __                         __  ###
### \  /  /\  |__) |  /\  |\ |  |     /  `  /\  |    |    | |\ | / _` ###
###  \/  /~~\ |  \ | /~~\ | \|  |     \__, /~~\ |___ |___ | | \| \__> ###
                                                                  
	#if not os.path.isdir('%s/tvc_de_novo' % intermediate_folder):
		#subprocess.call(['mkdir', '%s/tvc_de_novo' % intermediate_folder])
	#if not os.path.isdir('%s/tvc_only_hotspot' % intermediate_folder):
		#subprocess.call(['mkdir', '%s/tvc_only_hotspot' % intermediate_folder])
	
	## RUN 1 : DE NOVO
	#logging.info("\t - [%s] variantCaller <de novo> ..." % (time.strftime("%H:%M:%S")))
	#cmd = subprocess.Popen([
		#'python', '%s/variantCaller/bin/variant_caller_pipeline.py' % pipeline_folder,
		#'--input-bam',       bamfile,
		#'--output-dir',      '%s/tvc_de_novo' % intermediate_folder,
		#'--reference-fasta', reference,
		#'--region-bed',      bed_merged,
		#'--primer-trim-bed', target_bed,
		#'--parameters-file', param,
		#'--error-motifs',	'%s/variantCaller/share/TVC/sse/motifset.txt' % pipeline_folder,
		#'--postprocessed-bam','%s/tvc_de_novo/processed.bam' % intermediate_folder,
		#], 
		#stdout=open('%s/tvc_de_novo/variant_caller_pipeline.stdout.txt' % intermediate_folder,'w'),
		#stderr=open('%s/tvc_de_novo/variant_caller_pipeline.stderr.txt' % intermediate_folder,'w'))
	#cmd.communicate()

	#subprocess.call(['gzip','-d','-c','%s/tvc_de_novo/TSVC_variants.vcf.gz' % intermediate_folder],stdout=open('%s/tvc_de_novo/TSVC_variants.vcf' % intermediate_folder,'w'))

	## Generate Variant Table ('alleles.xls' file)
	#cmd = subprocess.Popen([
		#'python', '%s/variantCaller/bin/generate_variant_tables.py' % pipeline_folder,
		#'--input-vcf',		'%s/tvc_de_novo/TSVC_variants.vcf' % intermediate_folder,
		#'--region-bed',		target_bed,
		#'--output-xls',		'%s/tvc_de_novo/output.xls' % intermediate_folder,
		#'--alleles2-xls',	'%s/tvc_de_novo/alleles.xls' % intermediate_folder
		#], 
		#stdout=open('%s/tvc_de_novo/generate_variant_tables.stdout.txt' % intermediate_folder,'w'), 
		#stderr=open('%s/tvc_de_novo/generate_variant_tables.stderr.txt' % intermediate_folder,'w'))
	#cmd.communicate()
	
	## RUN 2 : ONLY HOTSPOT
	#logging.info("\t - [%s] variantCaller <only hotspot> ..." % (time.strftime("%H:%M:%S")))
	#if hotspot_vcf != '':
		#cmd = subprocess.Popen([
			#'python', '%s/variantCaller/bin/variant_caller_pipeline.py' % pipeline_folder,
			#'--input-bam',       bamfile,
			#'--output-dir',      '%s/tvc_only_hotspot' % intermediate_folder,
			#'--reference-fasta', reference,
			#'--region-bed',      bed_merged,
			#'--primer-trim-bed', target_bed,
			#'--parameters-file', param_hotspot_only,
			#'--error-motifs',	'%s/variantCaller/share/TVC/sse/motifset.txt' % pipeline_folder,
			#'--hotspot-vcf',     hotspot_vcf,
			#], 
			#stdout=open('%s/tvc_only_hotspot/variant_caller_pipeline.stdout.txt' % intermediate_folder,'w'),
			#stderr=open('%s/tvc_only_hotspot/variant_caller_pipeline.stderr.txt' % intermediate_folder,'w'))
		#cmd.communicate()	

		#subprocess.call(['gzip','-d','-c','%s/tvc_only_hotspot/TSVC_variants.vcf.gz' % intermediate_folder],stdout=open('%s/tvc_only_hotspot/TSVC_variants.vcf' % intermediate_folder,'w'))

		## Generate Variant Table ('alleles.xls' file)
		#cmd = subprocess.Popen([
			#'python', '%s/variantCaller/bin/generate_variant_tables.py' % pipeline_folder,
			#'--input-vcf',		'%s/tvc_only_hotspot/TSVC_variants.vcf' % intermediate_folder,
			#'--region-bed',		target_bed,
			#'--hotspots',
			#'--output-xls',		'%s/tvc_only_hotspot/output.xls' % intermediate_folder,
			#'--alleles2-xls',	'%s/tvc_only_hotspot/alleles.xls' % intermediate_folder
			#], 
			#stdout=open('%s/tvc_only_hotspot/generate_variant_tables.stdout.txt' % intermediate_folder,'w'), 
			#stderr=open('%s/tvc_only_hotspot/generate_variant_tables.stderr.txt' % intermediate_folder,'w'))
		#cmd.communicate()

####         __   ___  __  ___     __   __  ### 
#### | |\ | /__` |__  |__)  |     |  \ |__) ### 
#### | | \| .__/ |___ |  \  |     |__/ |__) ### 
                                          
	## add variantmetrics for each variant (if analysis already done, delete old variantmetrics first)
	## add new entries for new variants in db
	#logging.info("\t - [%s] insert new variants and metrics into DB ..." % (time.strftime("%H:%M:%S")))
	#abl1 = 'no'
	#if 'ABL1_NM_005157.fasta' in reference:
		#abl1 = 'yes'
		
	#cmd = subprocess.Popen(['python','%s/variantBase/insert_db_variants.py' % pipeline_folder,
		#'--analysis', analysis_id,
		#'--variants', '%s/tvc_de_novo/alleles.xls' % intermediate_folder,
		#'--abl1', abl1
		#], 
		#stdout=open('%s/insert_db_variants.stdout.txt' % intermediate_folder,'w'), 
		#stderr=open('%s/insert_db_variants.stderr.txt' % intermediate_folder,'w'))
	#cmd.communicate()
	#check_stderr('%s/insert_db_variants.stdout.txt'%intermediate_folder,indent=2)
	#check_stderr('%s/insert_db_variants.stderr.txt'%intermediate_folder,indent=2)
		
	####                 __  ___      ___  ___          ___                    __              ___  __  ###
	####  /\  |\ | |\ | /  \  |   /\   |  |__     |\ | |__  |  |    \  /  /\  |__) |  /\  |\ |  |  /__` ###
	#### /~~\ | \| | \| \__/  |  /~~\  |  |___    | \| |___ |/\|     \/  /~~\ |  \ | /~~\ | \|  |  .__/ ###

	#logging.info("\t - [%s] annotate new DB variants ..." % (time.strftime("%H:%M:%S")))
	
	#logging.info("\t\t - [%s] HGVS check ..." % (time.strftime("%H:%M:%S")))
	#cmd = subprocess.Popen(['python','%s/variantAnnotation/annotate_variantbase.step0.py' % pipeline_folder,'--new'],stdout=open('%s/annotate_variantbase.step0.stdout.txt' % intermediate_folder,'w'), stderr=open('%s/annotate_variantbase.step0.stderr.txt' % intermediate_folder,'w'))
	#cmd.communicate()
	#check_stderr('%s/annotate_variantbase.step0.stderr.txt' % intermediate_folder,indent=3)
	
	#logging.info("\t\t - [%s] Annovar and VEP ..." % (time.strftime("%H:%M:%S")))
	#cmd = subprocess.Popen(['python','%s/variantAnnotation/annotate_variantbase.step1.py' % pipeline_folder,'--new','--output-folder',intermediate_folder], stdout=open('%s/annotate_variantbase.step1.stdout.txt' % intermediate_folder,'w'), stderr=open('%s/annotate_variantbase.step1.stderr.txt' % intermediate_folder,'w'))
	#cmd.communicate()
	#check_stderr('%s/annotate_variantbase.step1.stderr.txt' % intermediate_folder,indent=3)

	#logging.info("\t\t - [%s] Merging annotations and updating DB ..." % (time.strftime("%H:%M:%S")))
	#if os.path.isfile('%s/annovar/annovar_input.tsv' % intermediate_folder):
		#cmd = subprocess.Popen(['python','%s/variantAnnotation/annotate_variantbase.step2.py' % pipeline_folder,'--annovar-results', '%s/annovar/annovar.hg19_multianno.txt' % intermediate_folder,'--vep-results', '%s/vep/vep_output.tsv' % intermediate_folder], stdout=open('%s/annotate_variantbase.step2.stdout.txt' % intermediate_folder,'w'),stderr=open('%s/annotate_variantbase.step2.stderr.txt' % intermediate_folder,'w'))
		#cmd.communicate()
		#check_stderr('%s/annotate_variantbase.step2.stderr.txt' % intermediate_folder,indent=3)
		
###  ___                   __   ___  __   __   __  ___ ###
### |__  | |\ |  /\  |    |__) |__  |__) /  \ |__)  |  ###
### |    | | \| /~~\ |___ |  \ |___ |    \__/ |  \  |  ###

	if not os.path.isdir('%s/finalReport' % intermediate_folder):
		subprocess.call(['mkdir', '%s/finalReport' % intermediate_folder])
	finalReport_path = '%s/%s_%s.finalReport.xlsx' % (sample_folder,sample,barcode)
	
	xmin = '300' # for coverage analysis sheet in finalreport
	if 'cDNA' in param_hotspot_only:
		xmin = '2000'
		
	logging.info("\t - [%s] finalReport ..." % (time.strftime("%H:%M:%S")))
	cmd = subprocess.Popen([
		'python','%s/finalReport/finalReport.py' % pipeline_folder,
		'--analysis', analysis_id,
		'--xmin', xmin,
	], 
	stdout=open('%s/finalReport/finalReport.stdout.txt' % intermediate_folder,'w'), 
	stderr=open('%s/finalReport/finalReport.stderr.txt' % intermediate_folder,'w'))
	cmd.communicate()
	check_stderr('%s/finalReport/finalReport.stdout.txt' % intermediate_folder,indent=2)
	check_stderr('%s/finalReport/finalReport.stderr.txt' % intermediate_folder,indent=2)

	#################################################
	
	if options.run:
		if 'HD802-HD748' in sample.upper():
			subprocess.call(['python','%s/scripts/HD802-HD748_check.py' % pipeline_folder,'%s/annovar/format_diag.tsv' % intermediate_folder,sample,run_name])
		elif 'ACROMETRIX' in sample.upper():
			subprocess.call(['python','%s/scripts/Acrometrix_check.py' % pipeline_folder,finalReport_path,sample,run_name])
		elif 'HORIZON' in sample.upper():
			subprocess.call(['python','%s/scripts/Horizon_check.py' % pipeline_folder,finalReport_path,sample,run_name])
		elif 'BARBI' in sample.upper():
			subprocess.call(['python','%s/scripts/Temoin_TP53_check.py' % pipeline_folder,finalReport_path,sample,run_name])
		for c in ['BAF-','BAF5','P190']:
			if c in sample.upper():
				subprocess.call(['python','%s/scripts/Temoins_ABL1_check.py' % pipeline_folder,finalReport_path,sample,run_name])
		if run_type == 'SBT':
			subprocess.call(['python','%s/scripts/collect_variants_MET_intron_13-14.py' % pipeline_folder,finalReport_path,sample,run_name])
		for control_name in control_names:
			if control_name in sample.upper():
				checkconta_bamlist.append(bamfile)

	subprocess.call(['cp',target_bed,intermediate_folder])
	subprocess.call(['cp',param,intermediate_folder])
	if os.path.isfile('%s/tvc_de_novo/processed.bam' % intermediate_folder):
		subprocess.call(['mv','%s/tvc_de_novo/processed.bam' % intermediate_folder,'%s/%s_%s.processed.bam' % (sample_folder,sample,barcode)])
	if os.path.isfile('%s/tvc_de_novo/processed.bam.bai' % intermediate_folder):
		subprocess.call(['mv','%s/tvc_de_novo/processed.bam.bai' % intermediate_folder,'%s/%s_%s.processed.bam.bai' % (sample_folder,sample,barcode)])
	
	shutil.make_archive(intermediate_folder,'zip',intermediate_folder)
	shutil.rmtree(intermediate_folder)
	
if options.run and os.path.isfile('%s/barcodes.json' % run_folder):
	
	### VBS scripts for printing
	logging.info(" [%s] Making VBS for fast printing..." % time.strftime("%H:%M:%S"))
	subprocess.call(['python','%s/finalReport/make_vbs_print.py' % pipeline_folder,run_folder,'%s/barcodes.json' % run_folder])
	
###  __        ___  __           __   __       ___                        ___    __       ###
### /  ` |__| |__  /  ` |__/    /  ` /  \ |\ |  |   /\   |\/| | |\ |  /\   |  | /  \ |\ | ###
### \__, |  | |___ \__, |  \    \__, \__/ | \|  |  /~~\  |  | | | \| /~~\  |  | \__/ | \| ###

	checkconta_folder = '%s/_checkContamination' % run_folder
	if not os.path.isdir(checkconta_folder):
		subprocess.call(['mkdir', checkconta_folder])
		
	for controlbam in checkconta_bamlist:
		logging.info(" [%s] checkContamination <%s> ..." % (time.strftime("%H:%M:%S"),controlbam.split('/')[-1].split('_IonXpress')[0]))
		cmd = subprocess.Popen(['python','%s/checkContamination/checkContamination.py' % pipeline_folder,'--bam', controlbam,'--read-len', str(checkconta_read_len)],stdout=open('%s/checkContamination.stdout.txt' % checkconta_folder,'w'),stderr=open('%s/checkContamination.stderr.txt' % checkconta_folder,'w'))
		cmd.communicate()
		check_stderr('%s/checkContamination.stderr.txt' % checkconta_folder,indent=1)
	
####  __        ___  __                  ___ ###
#### /  ` |__| |__  /  ` |__/  |\/| |  |  |  ###
#### \__, |  | |___ \__, |  \  |  | \__/  |  ###

	#checkMut_folder = '%s/_checkMut' % run_folder
	#if not os.path.isdir(checkMut_folder):
		#subprocess.call(['mkdir',checkMut_folder])
	#logging.info(" [%s] checkMut (routine variants) ..." % time.strftime("%H:%M:%S"))
	#subprocess.call(['python','%s/checkMut/routine_checkMut.py' % pipeline_folder,run_folder]) # ,stdout=open('%s/_routine_checkMut.stdout.txt' % checkMut_folder,'w'),stderr=open('%s/_routine_checkMut.stderr.txt' % checkMut_folder,'w')
	##cmd.communicate()
	##check_stderr('%s/_routine_checkMut.stderr.txt' % checkMut_folder,indent=1)
	
####  __        __   __  ###
#### |  \ |  | |__) / _` ###
#### |__/ \__/ |    \__> ###

	#dupG_folder = '%s/_dupG' % run_folder
	#if not os.path.isdir(dupG_folder):
		#subprocess.call(['mkdir',dupG_folder])
	#logging.info(" [%s] dupG.py ..." % time.strftime("%H:%M:%S"))
	#cmd = subprocess.Popen(['python', '%s/dupG/dupG.py' % pipeline_folder,'--run-folder', run_folder],stdout=open('%s/dupG.stdout.txt' % dupG_folder,'w'),stderr=open('%s/dupG.stderr.txt' % dupG_folder,'w'))	
	#cmd.communicate()
	#check_stderr('%s/dupG.stderr.txt' % dupG_folder,indent=1)
