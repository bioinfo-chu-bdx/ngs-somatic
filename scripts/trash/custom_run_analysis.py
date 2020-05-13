#!/usr/bin/python
import sys
import os
import time
import glob
import json
import subprocess
import shutil
from optparse import OptionParser

### GATHERING PARAMETERS ############################################################

parser = OptionParser()
parser.add_option('-b', '--bam',		help="Input bam file for SINGLE BAM ANALYSIS (must contain flow signal from TSS)", dest='bam')
parser.add_option('-r', '--full-run',	help="Run folder path  for FULL RUN ANALYSIS", dest='full_run') 
parser.add_option('-t', '--run-type',	help="Run type : SBT, LAM, TP53-HEMATO... (OPTIONAL, NOT NEEDED if barcodes.json is present in run folder)", dest='run_type') 
(options, args) = parser.parse_args()

if options.bam and options.full_run:
	sys.stderr.write("[run_analysis.py] Error: <--bam> and <--full-run> are not compatibles\n")
	sys.exit()
if options.full_run:
	bamlist = glob.glob(options.full_run+'/*/*.bam')
	bamlist = [item for item in bamlist if not 'processed' in item]
elif options.bam:
	bamlist = [options.bam]
else:
	sys.stderr.write("[run_analysis.py] Error: no <--bam> or <--full-run> specified\n")
	sys.exit()
	
with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)
fp_file = global_param['fp_file']
control_names = ['H2O','H20','NTC'] # liste des noms possibles pour les temoins negatifs
checkconta_bamlist = []
sampleask = False

### RUN ANALYSIS ######################################################################

bam_data = {}
	
for bamfile in sorted(bamlist) :
	sample = bamfile.split('/')[-1].split('_IonXpress')[0]
	barcode = 'IonXpress_' + bamfile.split('IonXpress_')[-1].split('.bam')[0]
	sample_folder = os.path.dirname(bamfile)
	run_folder = os.path.dirname(sample_folder)
	
	### RUN TYPE ###	
	run_type = None
	# if run-type manually specified
	if options.run_type:
		run_type = options.run_type
		if run_type not in global_param['run_type'].keys():
			sys.stderr.write("[run_analysis.py] Error: run type <%s> not found in global_parameters.json file\n" % run_type)
			sys.exit()
	# default : via barcodes.json in run folder
	elif os.path.isfile(run_folder+'/barcodes.json'):
		with open(run_folder+'/barcodes.json', 'r') as g:
			barcodes_json = json.load(g)
			bed_name = barcodes_json[barcode]['target_region_filepath'].split('/unmerged/detail/')[-1]
			for _run_type in global_param['run_type']:
				if global_param['run_type'][_run_type]['target_bed'].split('/')[-1] == bed_name:
					run_type = _run_type
					break
					
	### RUN TYPE PARAMETERS ###	
	if not run_type:
		print "\t -- Error : run type not found. Sample %s will not be processed." % sample
		continue
	reference = global_param['run_type'][run_type]['reference']
	target_bed = global_param['run_type'][run_type]['target_bed']
	bed_merged = global_param['run_type'][run_type]['merged_bed']
	param = global_param['run_type'][run_type]['vc_parameters']
	param_hotspot_only = global_param['run_type'][run_type]['vc_parameters_hotspot_only']
	dna_number = sample.split('-')[-1]
	if 'PL' in dna_number : # cDNA sample
		if 'vc_parameters_hotspot_cdna' in global_param['run_type'][run_type]:
			print "(sample %s is cDNA)" % sample
			param_hotspot_only = global_param['run_type'][run_type]['vc_parameters_hotspot_cdna']
	hotspot_vcf = global_param['run_type'][run_type]['hotspot_vcf']

	##### TO UNZIP ALL intermediates_files.zip (RERUN) ###
	#try:
		#subprocess.call(['unzip', '-o', intermediate_folder+'.zip','-d', intermediate_folder])
	#except:
		#print "* no intermediates_files.zip ton unzip"
	######################################################

	# dossier results
	intermediate_folder = sample_folder + '/intermediate_files'
	if not os.path.isdir(intermediate_folder):
		subprocess.call(['mkdir', intermediate_folder])
	
	bam_data[bamfile] = {'barcode':barcode,'sample':sample,'run_type':run_type,'reference':reference,'sample_folder':sample_folder,'intermediate_folder':intermediate_folder,'target_bed':target_bed,'bed_merged':bed_merged,'param':param,'param_hotspot_only':param_hotspot_only,'hotspot_vcf':hotspot_vcf}
		
###  __   __        ___  __        __   ___                             __     __  ###
### /  ` /  \ \  / |__  |__)  /\  / _` |__      /\  |\ |  /\  |    \ / /__` | /__` ###
### \__, \__/  \/  |___ |  \ /~~\ \__> |___    /~~\ | \| /~~\ |___  |  .__/ | .__/ ###
                                                                               
print " [%s] Coverage Analysis ..." % (time.strftime("%H:%M:%S"))
for bamfile in bam_data:
	coverage_folder = '%s/coverage' % bam_data[bamfile]['intermediate_folder']
	if not os.path.isdir(coverage_folder):
		subprocess.call(['mkdir', coverage_folder])

	cmd = subprocess.Popen([
		'bash', '/DATA/work/coverageAnalysis/run_coverage_analysis.sh',
		'-ag',
		'-D',coverage_folder,
		'-B',bam_data[bamfile]['target_bed'],
		bam_data[bamfile]['reference'],
		bamfile
		], stdout=open(coverage_folder+'/run_coverage_analysis.stdout.txt','w'), stderr=open(coverage_folder+'/run_coverage_analysis.stderr.txt','w'))	
		#		'-L', 'hg19',
	cmd.communicate()
	
###  __        __  ___  __   __        ___  __        __   ___ ###
### |__) |    /  \  |  /  ` /  \ \  / |__  |__)  /\  / _` |__  ###
### |    |___ \__/  |  \__, \__/  \/  |___ |  \ /~~\ \__> |___ ###

if options.full_run and os.path.isfile(run_folder+'/barcodes.json'):
	print " [%s] plotCoverage.py ..." % time.strftime("%H:%M:%S")
	cmd = subprocess.Popen(['python', '/DATA/work/plotCoverage/plotCoverage.py','--run-folder', run_folder], stdout=subprocess.PIPE)	
	cmd.communicate()

###  __            ###
### /  ` |\ | \  / ###
### \__, | \|  \/  ###

if options.full_run:
	if not os.path.isdir(run_folder+'/_CNA'):
		subprocess.call(['mkdir', run_folder+'/_CNA'])
	print " [%s] CNV Analysis ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen(['python','/DATA/work/CNV/run_cna.py', run_folder], stdout=open(run_folder+'/_CNA/run_cna.stdout.txt','w'), stderr=open(run_folder+'/_CNA/run_cna.stderr.txt','w'))


### RESTE DE L'ANALYSE PATIENT PAR PATIENT ###	
for bamfile in sorted(bamlist) :
	print " [%s] Processing %s ..." % (time.strftime("%H:%M:%S"),bamfile.split('/')[-1])
	
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
	
	print "\t * panel %s" % run_type	
	barcode = bam_data[bamfile]['barcode']
	sample = bam_data[bamfile]['sample']
	run_type = bam_data[bamfile]['run_type']
	reference = bam_data[bamfile]['reference']
	sample_folder = bam_data[bamfile]['sample_folder']
	intermediate_folder = bam_data[bamfile]['intermediate_folder']
	target_bed = bam_data[bamfile]['target_bed']
	bed_merged = bam_data[bamfile]['bed_merged']
	param = bam_data[bamfile]['param']
	param_hotspot_only = bam_data[bamfile]['param_hotspot_only']
	hotspot_vcf = bam_data[bamfile]['hotspot_vcf']
	
###            __              ___     __                         __  ###
### \  /  /\  |__) |  /\  |\ |  |     /  `  /\  |    |    | |\ | / _` ###
###  \/  /~~\ |  \ | /~~\ | \|  |     \__, /~~\ |___ |___ | | \| \__> ###
                                                                  
	if not os.path.isdir(intermediate_folder+'/tvc_de_novo'):
		subprocess.call(['mkdir', intermediate_folder+'/tvc_de_novo'])
	if not os.path.isdir(intermediate_folder+'/tvc_only_hotspot'):
		subprocess.call(['mkdir', intermediate_folder+'/tvc_only_hotspot'])
	
	#print 'python '+'/DATA/work/variantCaller/bin/variant_caller_pipeline.py '+'--input-bam '+bamfile+' --output-dir '+intermediate_folder+'/tvc_de_novo '+'--reference-fasta '+reference+' --region-bed '+bed_merged+' --primer-trim-bed '+target_bed+' --parameters-file '+param+' --error-motifs '+'/DATA/work/variantCaller/share/TVC/sse/motifset.txt '+'--postprocessed-bam '+intermediate_folder+'/tvc_de_novo'+'/processed.bam'
		
	# RUN 1 : DE NOVO
	print "\t - [%s] variant_caller_pipeline.py <de novo> ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen([
		'python', '/DATA/work/variantCaller/bin/variant_caller_pipeline.py',
		'--input-bam',       bamfile,
		'--output-dir',      intermediate_folder+'/tvc_de_novo',
		'--reference-fasta', reference,
		'--region-bed',      bed_merged,
		'--primer-trim-bed', target_bed,
		'--parameters-file', param,
		'--error-motifs',	'/DATA/work/variantCaller/share/TVC/sse/motifset.txt',
		'--postprocessed-bam',intermediate_folder+'/tvc_de_novo'+'/processed.bam',
		], stdout=open(intermediate_folder+'/tvc_de_novo'+'/variant_caller_pipeline.stdout.txt','w'),stderr=open(intermediate_folder+'/tvc_de_novo'+'/variant_caller_pipeline.stderr.txt','w'))
	cmd.communicate()

	subprocess.call(['gzip','-d','-c',intermediate_folder+'/tvc_de_novo'+'/TSVC_variants.vcf.gz'],stdout=open(intermediate_folder+'/tvc_de_novo'+'/TSVC_variants.vcf','w'))

	print "\t - [%s] generate_variant_tables.py ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen([
		'python', '/DATA/work/variantCaller/bin/generate_variant_tables.py',
		'--input-vcf',		intermediate_folder+'/tvc_de_novo'+'/TSVC_variants.vcf',
		'--region-bed',		target_bed,
		'--output-xls',		intermediate_folder+'/tvc_de_novo'+'/output.xls',
		'--alleles2-xls',	intermediate_folder+'/tvc_de_novo'+'/alleles.xls'
		], stdout=subprocess.PIPE)
	cmd.communicate()
	
	# RUN 2 : ONLY HOTSPOT
	print "\t - [%s] variant_caller_pipeline.py <only hotspot> ..." % (time.strftime("%H:%M:%S"))
	if hotspot_vcf != '':
		cmd = subprocess.Popen([
			'python', '/DATA/work/variantCaller/bin/variant_caller_pipeline.py',
			'--input-bam',       bamfile,
			'--output-dir',      intermediate_folder+'/tvc_only_hotspot',
			'--reference-fasta', reference,
			'--region-bed',      bed_merged,
			'--primer-trim-bed', target_bed,
			'--parameters-file', param_hotspot_only,
			'--error-motifs',	'/DATA/work/variantCaller/share/TVC/sse/motifset.txt',
			'--hotspot-vcf',     hotspot_vcf,
			], stdout=open(intermediate_folder+'/tvc_only_hotspot'+'/variant_caller_pipeline.stdout.txt','w'),stderr=open(intermediate_folder+'/tvc_only_hotspot'+'/variant_caller_pipeline.stderr.txt','w'))
		cmd.communicate()

		subprocess.call(['gzip','-d','-c',intermediate_folder+'/tvc_only_hotspot'+'/TSVC_variants.vcf.gz'],stdout=open(intermediate_folder+'/tvc_only_hotspot'+'/TSVC_variants.vcf','w'))

		print "\t - [%s] generate_variant_tables.py ..." % (time.strftime("%H:%M:%S"))
		cmd = subprocess.Popen([
			'python', '/DATA/work/variantCaller/bin/generate_variant_tables.py',
			'--input-vcf',		intermediate_folder+'/tvc_only_hotspot'+'/TSVC_variants.vcf',
			'--region-bed',		target_bed,
			'--hotspots',
			'--output-xls',		intermediate_folder+'/tvc_only_hotspot'+'/output.xls',
			'--alleles2-xls',	intermediate_folder+'/tvc_only_hotspot'+'/alleles.xls'
			], stdout=subprocess.PIPE)
		cmd.communicate()

###                 __  ___      ___    __       ###
###  /\  |\ | |\ | /  \  |   /\   |  | /  \ |\ | ###
### /~~\ | \| | \| \__/  |  /~~\  |  | \__/ | \| ###

#####                 __             __  #####
#####  /\  |\ | |\ | /  \ \  /  /\  |__) #####
##### /~~\ | \| | \| \__/  \/  /~~\ |  \ #####
	if not os.path.isdir(intermediate_folder+'/annovar'):
		subprocess.call(['mkdir', intermediate_folder+'/annovar'])
			
	print "\t - [%s] generate_annovar_input.py ..." % (time.strftime("%H:%M:%S"))
	if 'ABL1_NM_005157.fasta' in reference:
		cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/custom_generate_annovar_input.py',
		'--input', intermediate_folder+'/tvc_de_novo/alleles.xls',
		'--output', intermediate_folder+'/annovar/annovar_input.tsv',
		'--abl1'
		], stdout=open(intermediate_folder+'/annovar/generate_annovar_input.stdout.txt','w'), stderr=open(intermediate_folder+'/annovar/generate_annovar_input.stderr.txt','w'))
		cmd.communicate()
	else:
		cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/custom_generate_annovar_input.py',
		'--input', intermediate_folder+'/tvc_de_novo/alleles.xls',
		'--output', intermediate_folder+'/annovar/annovar_input.tsv'
		], stdout=open(intermediate_folder+'/annovar/generate_annovar_input.stdout.txt','w'), stderr=open(intermediate_folder+'/annovar/generate_annovar_input.stderr.txt','w'))
		cmd.communicate()

	print "\t - [%s] table_annovar.pl ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen([
		'perl', '/DATA/work/variantAnnotation/annovar/table_annovar.pl', 
		intermediate_folder+'/annovar/annovar_input.tsv', '/DATA/work/variantAnnotation/annovar/humandb/',
		'-buildver', 'hg19',
		'-out', intermediate_folder+'/annovar/variants.annotated',
		'-protocol', 'refGeneWithVer,cosmic89,avsnp150,intervar_20180118,clinvar_20190305,nci60,esp6500siv2_all,1000g2015aug_all,gnomad211_genome,exac03,dbnsfp33a',
		'-operation', 'g,f,f,f,f,f,f,f,f,f,f',
		'-argument', '--hgvs,,,,,,,,,,',
		'-nastring', '.',
		'-polish',
		'-xref', '/DATA/work/variantAnnotation/annovar/example/gene_xref.txt',
		], stdout=open(intermediate_folder+'/annovar/table_annovar.stdout.txt','w'), stderr=open(intermediate_folder+'/annovar/table_annovar.stderr.txt','w'))
	cmd.communicate()
#####       ___  __  #####
##### \  / |__  |__) #####
#####  \/  |___ |    #####
	if not os.path.isdir(intermediate_folder+'/vep'):
		subprocess.call(['mkdir', intermediate_folder+'/vep'])
	
	if 'ABL1_NM_005157.fasta' in reference:
		print "\t - [%s] generate_vep_input.py ..." % (time.strftime("%H:%M:%S"))
		cmd = subprocess.Popen([
			'python', '/DATA/work/variantAnnotation/custom_generate_vep_input.py',
			'--input', intermediate_folder+'/tvc_de_novo/alleles.xls',
			'--output', intermediate_folder+'/vep/vep_input.tsv',
			'--abl1'
			], stdout=open(intermediate_folder+'/vep/generate_vep_input.stdout.txt','w'), stderr=open(intermediate_folder+'/vep/generate_vep_input.stderr.txt','w'))
		cmd.communicate()
	else:
		print "\t - [%s] generate_vep_input.py ..." % (time.strftime("%H:%M:%S"))
		cmd = subprocess.Popen([
			'python', '/DATA/work/variantAnnotation/custom_generate_vep_input.py',
			'--input', intermediate_folder+'/tvc_de_novo/alleles.xls',
			'--output', intermediate_folder+'/vep/vep_input.tsv'
			], stdout=open(intermediate_folder+'/vep/generate_vep_input.stdout.txt','w'), stderr=open(intermediate_folder+'/vep/generate_vep_input.stderr.txt','w'))
		cmd.communicate()

	print "\t - [%s] run_vep ..." % (time.strftime("%H:%M:%S"))
	#cmd = subprocess.Popen(['python', '/DATA/work/variantAnnotation/VEP/run_vep.py', intermediate_folder+'/vep/vep_input.tsv'], stdout=open(intermediate_folder+'/vep/vep.stdout.txt','w'), stderr=open(intermediate_folder+'/vep/vep.stderr.txt','w'))
	cmd = subprocess.Popen([
		'perl', '/DATA/work/variantAnnotation/VEP/ensembl-vep/vep',
		'--offline',
		'--dir_cache', '/DATA/work/variantAnnotation/VEP/cache/',
		'--force_overwrite',
		'--refseq',
		'--numbers',
		'--hgvs',
		'--hgvsg',
		'--no_escape',
		'--check_existing', # retire des "existing variants" les variants connus qui ne correspondent pas aux alleles mutes (par defaut utile seulement les coordonnees)
		# For some data sources (COSMIC, HGMD), Ensembl is not licensed to redistribute allele-specific data, so VEP will report the existence of co-located variants with unknown alleles without carrying out allele matching.
		'--variant_class',
		'--sift', 'p',
		'--polyphen', 'p',
		'--af_1kg',
		'--af',
		'--freq_pop', '1KG_ALL',
		'--freq_pop', 'gnomAD',
		'--af_esp',
		'--symbol',
		'--tab',
		'--pubmed',
		'--fasta', '/DATA/work/variantAnnotation/VEP/cache/homo_sapiens/96_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz',
		'--exclude_predicted', # RETIRE LES TRANSCRITS PREDITS "XM" ou "XR"
		'--no_stats',
		'--verbose',
		'-i', intermediate_folder+'/vep/vep_input.tsv',
		'-o', intermediate_folder+'/vep/vep_output.tsv'
		], stdout=open(intermediate_folder+'/vep/run_vep.stdout.txt','w'), stderr=open(intermediate_folder+'/vep/run_vep.stderr.txt','w'))
	
	cmd.communicate()
#####  ___  __   __             ___     __          __  #####
##### |__  /  \ |__)  |\/|  /\   |     |  \ |  /\  / _` #####
##### |    \__/ |  \  |  | /~~\  |     |__/ | /~~\ \__> #####                                         
	print "\t - [%s] format_diag.py ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/custom_format_diag.py',
		'--annovar-input', intermediate_folder+'/annovar/variants.annotated.hg19_multianno.txt',
		'--vep-input', intermediate_folder+'/vep/vep_output.tsv',
		'--bed', target_bed
		], stdout=open(intermediate_folder+'/annovar/format_diag.stdout.txt','w'), stderr=open(intermediate_folder+'/annovar/format_diag.stderr.txt','w'))
	cmd.communicate()

###  ___                   __   ___  __   __   __  ___ ###
### |__  | |\ |  /\  |    |__) |__  |__) /  \ |__)  |  ###
### |    | | \| /~~\ |___ |  \ |___ |    \__/ |  \  |  ###

	if not os.path.isdir(intermediate_folder+'/finalReport'):
		subprocess.call(['mkdir', intermediate_folder+'/finalReport'])
	finalreport_name = 'test_new_'+sample+'_'+barcode+'.finalReport.xlsx'
	finalReport_path = intermediate_folder+'/finalReport/'+finalreport_name
	
	xmin = '300' # for coverage analysis sheet in finalreport
	if 'cDNA' in param_hotspot_only:
		xmin = '2000'
		
	print "\t - [%s] finalReport.py ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen([
		'python','/DATA/work/finalReport/custom_finalReport.py',
		'--bam', bamfile,
		'--run-type', run_type,
		'--output', finalReport_path,
		'--xmin', xmin,
	], stdout=open(intermediate_folder+'/finalReport/finalReport.stdout.txt','w'), stderr=open(intermediate_folder+'/finalReport/finalReport.stderr.txt','w'))
	cmd.communicate()

	#################################################
	
	if run_type == 'Lymphome_B' or run_type == 'Lymphome_T':
		subprocess.call(['python','/DATA/work/scripts/annot_lymphT_temoincnv_efs.py','--bam',bamfile])

	if options.full_run:
		run_name = run_folder.split('/')[-1]
		if 'HD802-HD748' in sample.upper(): # HD802-HD748 obsolete?
			subprocess.call(['python','/DATA/work/scripts/HD802-HD748_check.py',intermediate_folder+'/annovar/format_diag.tsv',sample,run_name])
		if 'ACROMETRIX' in sample.upper():
			subprocess.call(['python','/DATA/work/scripts/Acrometrix_check.py',finalReport_path,sample,run_name])
		if 'HORIZON' in sample.upper():
			subprocess.call(['python','/DATA/work/scripts/Horizon_check.py',finalReport_path,sample,run_name])
		if 'BAF-' in sample.upper():
			subprocess.call(['python','/DATA/work/scripts/Temoin_ABL1_check.py',finalReport_path,sample,run_name])
		if run_type == 'SBT':
			subprocess.call(['python','/DATA/work/scripts/collect_variants_MET_intron_13-14.py',finalReport_path,sample,run_name])
		for control_name in control_names:
			if control_name in sample.upper():
				checkconta_bamlist.append(bamfile)

	subprocess.call(['cp',target_bed,intermediate_folder])
	subprocess.call(['cp',param,intermediate_folder])
	subprocess.call(['mv',finalReport_path,sample_folder])
	subprocess.call(['mv',intermediate_folder+'/tvc_de_novo/processed.bam','%s/%s_%s.processed.bam' % (sample_folder,sample,barcode)])
	subprocess.call(['mv',intermediate_folder+'/tvc_de_novo/processed.bam.bai','%s/%s_%s.processed.bam.bai' % (sample_folder,sample,barcode)])
	subprocess.call(['mv',intermediate_folder+'/test_new_ALAMUT_load_variants.vbs',sample_folder])
	subprocess.call(['mv',intermediate_folder+'/test_new_ALAMUT_load_processed_bam.vbs',sample_folder])
	subprocess.call(['mv',intermediate_folder+'/test_new_PRINT_finalReport.vbs',sample_folder])
	
	shutil.make_archive(intermediate_folder,'zip',intermediate_folder)
	shutil.rmtree(intermediate_folder)

	print " [%s] ... Done" % (time.strftime("%H:%M:%S"))
	
	
if options.full_run and os.path.isfile(run_folder+'/barcodes.json'):
	
	### VBS scripts for printing
	print " [%s] Making VBS for fast printing..." % (time.strftime("%H:%M:%S"))
	subprocess.call(['python','/DATA/work/finalReport/custom_make_vbs_print.py',run_folder,run_folder+'/barcodes.json'])
	
###  __        ___  __           __   __       ___                        ___    __       ###
### /  ` |__| |__  /  ` |__/    /  ` /  \ |\ |  |   /\   |\/| | |\ |  /\   |  | /  \ |\ | ###
### \__, |  | |___ \__, |  \    \__, \__/ | \|  |  /~~\  |  | | | \| /~~\  |  | \__/ | \| ###
                                                                                      
	checkconta_folder = run_folder+'/_Check-contamination'
	if not os.path.isdir(checkconta_folder):
		subprocess.call(['mkdir', checkconta_folder])
		
	for controlbam in checkconta_bamlist:
		print " [%s] Check_contamination.py <%s> ..." % (time.strftime("%H:%M:%S"),controlbam.split('/')[-1].split('_IonXpress')[0])
		cmd = subprocess.Popen(['python','/DATA/work/checkContamination/custom_Check_contamination.py','--bam', controlbam,'--read-len', '100',], stdout=subprocess.PIPE)
		cmd.communicate()
	
###  __        ___  __                  ___ ###
### /  ` |__| |__  /  ` |__/  |\/| |  |  |  ###
### \__, |  | |___ \__, |  \  |  | \__/  |  ###
                                        
	#print " [%s] checkMut.py ..." % time.strftime("%H:%M:%S")
	#subprocess.call(['python','/DATA/work/scripts/routine_checkMut.py',options.full_run])
	
###  __        __   __  ###
### |  \ |  | |__) / _` ###
### |__/ \__/ |    \__> ###

	#print " [%s] dupG.py ..." % time.strftime("%H:%M:%S")
	#cmd = subprocess.Popen(['python', '/DATA/work/dupG/dupG.py','--run-folder', run_folder], stdout=subprocess.PIPE)	
	#cmd.communicate()

###            __              ___  __        __   ___ ###
### \  /  /\  |__) |  /\  |\ |  |  |__)  /\  /__` |__  ###
###  \/  /~~\ |  \ | /~~\ | \|  |  |__) /~~\ .__/ |___ ###
                                                   	
	#print " [%s] variantCollector.py ..." % time.strftime("%H:%M:%S")
	#cmd = subprocess.Popen(['python', '/DATA/work/variantBase/variantCollector.py','--run-folder', run_folder,'--sbt-wait'], preexec_fn=os.setpgrp)
	#cmd.communicate()
