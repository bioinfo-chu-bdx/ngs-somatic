#!/usr/bin/python
import sys
import os
import time
import glob
import json
import subprocess
import zipfile
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

	# dossier results
	intermediate_folder = sample_folder + '/intermediate_files'
	if not os.path.isdir(intermediate_folder):
		subprocess.call(['mkdir', intermediate_folder])
	
	bam_data[bamfile] = {'barcode':barcode,'sample':sample,'run_type':run_type,'reference':reference,'sample_folder':sample_folder,'intermediate_folder':intermediate_folder,'target_bed':target_bed,'bed_merged':bed_merged,'param':param,'param_hotspot_only':param_hotspot_only,'hotspot_vcf':hotspot_vcf}
		
#######################
## COVERAGE ANALYSIS ## 
#######################

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
	
	###################
	## VARIANTCALLER ##
	###################
	
	if not os.path.isdir(intermediate_folder+'/tvc_de_novo'):
		subprocess.call(['mkdir', intermediate_folder+'/tvc_de_novo'])
	if not os.path.isdir(intermediate_folder+'/tvc_only_hotspot'):
		subprocess.call(['mkdir', intermediate_folder+'/tvc_only_hotspot'])
	
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
			
	#######################
	## VARIANTANNOTATION ##
	#######################

	# ANNOVAR
	if not os.path.isdir(intermediate_folder+'/annovar'):
		subprocess.call(['mkdir', intermediate_folder+'/annovar'])
			
	print "\t - [%s] generate_annovar_input.py ..." % (time.strftime("%H:%M:%S"))
	if 'ABL1_NM_005157.fasta' in reference:
		cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/generate_annovar_input.py',intermediate_folder+'/tvc_de_novo/alleles.xls',intermediate_folder+'/annovar/annovar_input.tsv',global_param['abl1_cdna2genomic']], stdout=open(intermediate_folder+'/annovar/generate_annovar_input.stdout.txt','w'), stderr=open(intermediate_folder+'/annovar/generate_annovar_input.stderr.txt','w'))
		cmd.communicate()
	else:
		cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/generate_annovar_input.py',intermediate_folder+'/tvc_de_novo/alleles.xls',intermediate_folder+'/annovar/annovar_input.tsv'], stdout=open(intermediate_folder+'/annovar/generate_annovar_input.stdout.txt','w'), stderr=open(intermediate_folder+'/annovar/generate_annovar_input.stderr.txt','w'))
		cmd.communicate()

	print "\t - [%s] table_annovar.pl ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen([
		'perl', '/DATA/work/variantAnnotation/annovar/table_annovar.pl', intermediate_folder+'/annovar/annovar_input.tsv', '/DATA/work/variantAnnotation/annovar/humandb/',
		'-buildver', 'hg19',
		'-out', intermediate_folder+'/annovar/variants.annotated',
		'-protocol', 'refGene,cosmic86,avsnp150,intervar_20180118,clinvar_20180603,nci60,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_sas,dbnsfp33a',
		'-operation', 'g,f,f,f,f,f,f,f,f,f,f,f,f,f',
		'-argument', '--hgvs,,,,,,,,,,,,,',
		'-nastring', '.',
		'-polish',
		'-xref', '/DATA/work/variantAnnotation/annovar/example/gene_xref.txt',
		], stdout=open(intermediate_folder+'/annovar/table_annovar.stdout.txt','w'), stderr=open(intermediate_folder+'/annovar/table_annovar.stderr.txt','w'))
	cmd.communicate()

	## format diag
	print "\t - [%s] format_diag.py ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/format_diag.py',intermediate_folder+'/annovar/variants.annotated.hg19_multianno.txt',target_bed], stdout=open(intermediate_folder+'/annovar/format_diag.stdout.txt','w'), stderr=open(intermediate_folder+'/annovar/format_diag.stderr.txt','w'))
	cmd.communicate()
	
	## VEP (all  transcripts)
	if not os.path.isdir(intermediate_folder+'/vep'):
		subprocess.call(['mkdir', intermediate_folder+'/vep'])
	
	if 'ABL1_NM_005157.fasta' in reference:
		print "\t - [%s] vep_abl1.py ..." % (time.strftime("%H:%M:%S"))
		cmd = subprocess.Popen(['python', '/DATA/work/variantAnnotation/VEP/vep_abl1.py', intermediate_folder+'/tvc_de_novo/alleles.xls', intermediate_folder+'/vep/format_vep.tsv'], stdout=open(intermediate_folder+'/vep/vep.stdout.txt','w'), stderr=open(intermediate_folder+'/vep/vep.stderr.txt','w'))
		cmd.communicate()
	else:
		print "\t - [%s] generate_vep_input.py ..." % (time.strftime("%H:%M:%S"))
		cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/generate_vep_input.py',intermediate_folder+'/tvc_de_novo/alleles.xls',intermediate_folder+'/vep/vep_input.tsv'], stdout=open(intermediate_folder+'/vep/generate_vep_input.stdout.txt','w'), stderr=open(intermediate_folder+'/vep/generate_vep_input.stderr.txt','w'))
		cmd.communicate()

		print "\t - [%s] run_vep.py ..." % (time.strftime("%H:%M:%S"))
		cmd = subprocess.Popen(['python', '/DATA/work/variantAnnotation/VEP/run_vep.py', intermediate_folder+'/vep/vep_input.tsv'], stdout=open(intermediate_folder+'/vep/vep.stdout.txt','w'), stderr=open(intermediate_folder+'/vep/vep.stderr.txt','w'))
		cmd.communicate()

	#################
	## FINALREPORT ##
	#################
	
	if not os.path.isdir(intermediate_folder+'/finalReport'):
		subprocess.call(['mkdir', intermediate_folder+'/finalReport'])
	finalreport_name = sample+'_'+barcode+'.finalReport.xlsx'
	finalReport_path = intermediate_folder+'/finalReport/'+finalreport_name
	
	print "\t - [%s] finalReport.py ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen([
		'python','/DATA/work/finalReport/finalReport.py',
		'--bam', bamfile,
		'--run-type', run_type,
		'--output', finalReport_path,
	], stdout=open(intermediate_folder+'/finalReport/finalReport.stdout.txt','w'), stderr=open(intermediate_folder+'/finalReport/finalReport.stderr.txt','w'))
	cmd.communicate()

	################################################

	subprocess.call(['cp',target_bed,intermediate_folder])
	subprocess.call(['cp',param,intermediate_folder])
	subprocess.call(['mv',finalReport_path,sample_folder])
	subprocess.call(['mv',intermediate_folder+'/_ALAMUT_load_variants.vbs',sample_folder])
	subprocess.call(['mv',intermediate_folder+'/_ALAMUT_load_processed_bam.vbs',sample_folder])
	subprocess.call(['mv',intermediate_folder+'/_PRINT_finalReport.vbs',sample_folder])
	
	shutil.make_archive(intermediate_folder,'zip',intermediate_folder)
	shutil.rmtree(intermediate_folder)

	print " [%s] ... Done" % (time.strftime("%H:%M:%S"))
	

