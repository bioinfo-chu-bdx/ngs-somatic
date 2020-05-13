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
reference = global_param['ref']
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
	if run_type:
		target_bed = global_param['run_type'][run_type]['target_bed']
		bed_merged = global_param['run_type'][run_type]['merged_bed']
		param = global_param['run_type'][run_type]['vc_parameters']
		param_hotspot_only = global_param['run_type'][run_type]['vc_parameters_hotspot_only']
		hotspot_vcf = global_param['run_type'][run_type]['hotspot_vcf']
	else:
		print "\t -- Error : run type not found. Sample %s will not be processed." % sample
		continue
	ds_file = str(global_param.get('details_sensi_file').get(run_type)) # str() car si None cela devient "None"
	# dossier results
	#results_name = bamfile.split('/')[-1].replace(' ','@') # necessaire, variant_caller_pipeline.py plante si espace dans output_dir
	intermediate_folder = sample_folder + '/intermediate_files'
	if zipfile.is_zipfile(intermediate_folder+'.zip'):
		z = zipfile.ZipFile(intermediate_folder+'.zip','r')
		z.extractall(intermediate_folder)
		z.close()
	elif not os.path.isdir(intermediate_folder):
		subprocess.call(['mkdir', intermediate_folder])
	
	bam_data[bamfile] = {'barcode':barcode,'sample':sample,'run_type':run_type,'sample_folder':sample_folder,'intermediate_folder':intermediate_folder,'target_bed':target_bed,'bed_merged':bed_merged,'param':param,'param_hotspot_only':param_hotspot_only,'hotspot_vcf':hotspot_vcf,'ds_file':ds_file}
		
#######################
## COVERAGE ANALYSIS ## 
#######################

#print " [%s] Coverage Analysis ..." % (time.strftime("%H:%M:%S"))	
#for bamfile in bam_data:
	#coverage_folder = '%s/coverage' % bam_data[bamfile]['intermediate_folder']
	#if not os.path.isdir(coverage_folder):
		#subprocess.call(['mkdir', coverage_folder])
		
	#cmd = subprocess.Popen([
		#'bash', '/DATA/work/coverageAnalysis/run_coverage_analysis.sh',
		#'-ag',
		#'-L', 'hg19',
		#'-D',coverage_folder,
		#'-B',bam_data[bamfile]['target_bed'],
		#reference,
		#bamfile
		#], stdout=subprocess.PIPE)	
	#out, err = cmd.communicate()

##########
### CNV ##
##########

#if options.full_run:
	#print " [%s] CNV Analysis ..." % (time.strftime("%H:%M:%S"))
	#cmd = subprocess.call(['python','/DATA/work/CNV/run_cna.py', run_folder])
	
	
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
	
	print "\t -run_type =  %s" % run_type	
	barcode = bam_data[bamfile]['barcode']
	sample = bam_data[bamfile]['sample']
	run_type = bam_data[bamfile]['run_type']
	sample_folder = bam_data[bamfile]['sample_folder']
	intermediate_folder = bam_data[bamfile]['intermediate_folder']
	target_bed = bam_data[bamfile]['target_bed']
	bed_merged = bam_data[bamfile]['bed_merged']
	param = bam_data[bamfile]['param']
	param_hotspot_only = bam_data[bamfile]['param_hotspot_only']
	hotspot_vcf = bam_data[bamfile]['hotspot_vcf']
	ds_file = bam_data[bamfile]['ds_file']
	
	###################
	## VARIANTCALLER ##
	####################
	
	#if not os.path.isdir(intermediate_folder+'/tvc_de_novo'):
		#subprocess.call(['mkdir', intermediate_folder+'/tvc_de_novo'])
	#if not os.path.isdir(intermediate_folder+'/tvc_only_hotspot'):
		#subprocess.call(['mkdir', intermediate_folder+'/tvc_only_hotspot'])
	
	## RUN 1 : DE NOVO
	
	#cmd = subprocess.Popen([
		#'python', '/DATA/work/variantCaller/bin/variant_caller_pipeline.py',
		#'--input-bam',       bamfile,
		#'--output-dir',      intermediate_folder+'/tvc_de_novo',
		#'--reference-fasta', reference,
		#'--region-bed',      bed_merged,
		#'--primer-trim-bed', target_bed,
		#'--parameters-file', param,
		#'--error-motifs',	'/DATA/work/variantCaller/share/TVC/sse/motifset.txt',
		#'--postprocessed-bam',intermediate_folder+'/tvc_de_novo'+'/processed.bam',
		#], stdout=open(intermediate_folder+'/tvc_de_novo'+'/variant_caller_pipeline.stdout.txt','w'),stderr=open(intermediate_folder+'/tvc_de_novo'+'/variant_caller_pipeline.stderr.txt','w'))
	#out, err = cmd.communicate()
	#print "\t -variant_caller_pipeline.py <de novo> done, OUT: %s, ERR: %s"% (out, err)

	#subprocess.call(['gzip','-d','-c',intermediate_folder+'/tvc_de_novo'+'/TSVC_variants.vcf.gz'],stdout=open(intermediate_folder+'/tvc_de_novo'+'/TSVC_variants.vcf','w'))

	#cmd = subprocess.Popen([
		#'python', '/DATA/work/variantCaller/bin/generate_variant_tables.py',
		#'--input-vcf',		intermediate_folder+'/tvc_de_novo'+'/TSVC_variants.vcf',
		#'--region-bed',		target_bed,
		#'--output-xls',		intermediate_folder+'/tvc_de_novo'+'/output.xls',
		#'--alleles2-xls',	intermediate_folder+'/tvc_de_novo'+'/alleles.xls'
		#], stdout=subprocess.PIPE)
	#out, err = cmd.communicate()
	#print "\t -generate_variant_tables.py <de novo> done, OUT: %s, ERR: %s"%(out, err)	
	
	## RUN 2 : ONLY HOTSPOT
	
	#if hotspot_vcf != '':
		#cmd = subprocess.Popen([
			#'python', '/DATA/work/variantCaller/bin/variant_caller_pipeline.py',
			#'--input-bam',       bamfile,
			#'--output-dir',      intermediate_folder+'/tvc_only_hotspot',
			#'--reference-fasta', reference,
			#'--region-bed',      bed_merged,
			#'--primer-trim-bed', target_bed,
			#'--parameters-file', param_hotspot_only,
			#'--error-motifs',	'/DATA/work/variantCaller/share/TVC/sse/motifset.txt',
			#'--hotspot-vcf',     hotspot_vcf,
			#], stdout=open(intermediate_folder+'/tvc_only_hotspot'+'/variant_caller_pipeline.stdout.txt','w'),stderr=open(intermediate_folder+'/tvc_only_hotspot'+'/variant_caller_pipeline.stderr.txt','w'))
		#out, err = cmd.communicate()
		#print "\t -variant_caller_pipeline.py <only hotspot> done, OUT: %s, ERR: %s"% (out, err)

		#subprocess.call(['gzip','-d','-c',intermediate_folder+'/tvc_only_hotspot'+'/TSVC_variants.vcf.gz'],stdout=open(intermediate_folder+'/tvc_only_hotspot'+'/TSVC_variants.vcf','w'))

		#cmd = subprocess.Popen([
			#'python', '/DATA/work/variantCaller/bin/generate_variant_tables.py',
			#'--input-vcf',		intermediate_folder+'/tvc_only_hotspot'+'/TSVC_variants.vcf',
			#'--region-bed',		target_bed,
			#'--hotspots',
			#'--output-xls',		intermediate_folder+'/tvc_only_hotspot'+'/output.xls',
			#'--alleles2-xls',	intermediate_folder+'/tvc_only_hotspot'+'/alleles.xls'
			#], stdout=subprocess.PIPE)
		#out, err = cmd.communicate()
		#print "\t -generate_variant_tables.py <only hotspot> done, OUT: %s, ERR: %s"%(out, err)	
			
	########################
	### VARIANTANNOTATION ##
	########################

	#cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/generate_annovar_avinput.py',intermediate_folder+'/tvc_de_novo/alleles.xls'], stdout=open(intermediate_folder+'/tvc_de_novo/generate_annovar_avinput.stdout.txt','w'), stderr=open(intermediate_folder+'/tvc_de_novo/generate_annovar_avinput.stderr.txt','w'))
	#out, err = cmd.communicate()
	#print "\t -generate_annovar_avinput.py <de novo> done, OUT: %s, ERR: %s"%(out, err)
	#if hotspot_vcf != '':
		#cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/generate_annovar_avinput.py',intermediate_folder+'/tvc_only_hotspot/alleles.xls'], stdout=open(intermediate_folder+'/tvc_only_hotspot/generate_annovar_avinput.stdout.txt','w'), stderr=open(intermediate_folder+'/tvc_only_hotspot/generate_annovar_avinput.stderr.txt','w'))
		#out, err = cmd.communicate()
		#print "\t -generate_annovar_avinput.py <only hotspot> done, OUT: %s, ERR: %s"%(out, err)
	#cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/merge_annovar_avinput.py',intermediate_folder], stdout=open(intermediate_folder+'/merge_annovar_avinput.stdout.txt','w'), stderr=open(intermediate_folder+'/merge_annovar_avinput.stderr.txt','w'))
	#out, err = cmd.communicate()
	#print "\t -merge_annovar_avinput.py done, OUT: %s, ERR: %s"%(out, err)

	#cmd = subprocess.Popen([
		#'perl', '/DATA/work/variantAnnotation/annovar/table_annovar.pl', intermediate_folder+'/variants.merged.avinput', '/DATA/work/variantAnnotation/annovar/humandb/',
		#'-buildver', 'hg19',
		#'-out', intermediate_folder+'/variants.annotated',
		#'-protocol', 'refGene,cosmic82,avsnp147,intervar_20170202,clinvar_20160302,nci60,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_sas,dbnsfp33a',
		#'-operation', 'g,f,f,f,f,f,f,f,f,f,f,f,f,f',
		#'-argument', '--hgvs,,,,,,,,,,,,,',
		#'-nastring', '.',
		#'-polish',
		#'-xref', '/DATA/work/variantAnnotation/annovar/example/gene_xref.txt',
		#], stdout=open(intermediate_folder+'/table_annovar.stdout.txt','w'), stderr=open(intermediate_folder+'/table_annovar.stderr.txt','w'))
	#out, err = cmd.communicate()
	#print "\t -table_annovar.pl done, ERR: %s"% err

	### format diag
	#cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/format_diag.py',intermediate_folder+'/variants.annotated.hg19_multianno.txt',target_bed], stdout=subprocess.PIPE)
	#out, err = cmd.communicate()
	#print "\t -format_diag.py done, OUT: %s, ERR: %s"%(out, err)
	
	## OPTIONAL : format VEP with all possible transcripts
	cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/generate_vep_input.py',intermediate_folder+'/tvc_de_novo/alleles.xls',intermediate_folder+'/vep_input.tsv'], stdout=open(intermediate_folder+'/tvc_de_novo/generate_annovar_avinput.stdout.txt','w'), stderr=open(intermediate_folder+'/tvc_de_novo/generate_annovar_avinput.stderr.txt','w'))
	out, err = cmd.communicate()
	print "\t -generate_annovar_avinput.py <de novo> done, OUT: %s, ERR: %s"%(out, err)
	
	cmd = subprocess.Popen(['python', '/DATA/work/variantAnnotation/VEP/vep_local.py', intermediate_folder+'/vep_input.tsv'], stdout=open(intermediate_folder+'/vep.stdout.txt','w'), stderr=open(intermediate_folder+'/vep.stderr.txt','w'))
	out, err = cmd.communicate()
	print "\t -vep_local.py done, ERR: %s"% err
	
	#cmd = subprocess.Popen(['python', '/DATA/work/variantAnnotation/VEP/vep.py', intermediate_folder+'/variants.merged.avinput', intermediate_folder+'/format_vep.tsv'], stdout=open(intermediate_folder+'/vep.stdout.txt','w'), stderr=open(intermediate_folder+'/vep.stderr.txt','w'))
	#out, err = cmd.communicate()
	#print "\t -vep.py done, ERR: %s"% err

	#################
	## FINALREPORT ##
	#################
	if run_type == 'old_SBT_v5':
		run_type = 'SBT'
	
	finalreport_name = sample+'_'+barcode+'.finalReport.xlsx'
	finalReport_path = intermediate_folder+'/'+finalreport_name
	cmd = subprocess.Popen([
		'python','/DATA/work/finalReport/finalReport.py',
		'--input-tsv', intermediate_folder+'/format_diag.tsv',
		'--output', finalReport_path,
		'--sample',	sample,
		'--barcode', barcode,
		'--sample_folder', sample_folder,
		'--run_type', run_type,
		'--target-bed', target_bed,
		'--fp-file', fp_file,
		'--ds-file', ds_file,
	], stdout=open(intermediate_folder+'/finalReport.stdout.txt','w'))
	out, err = cmd.communicate()
	print "\t -finalReport.py done, ERR: %s"% err

	################################################
	# copy finalreport and vbs to bamfile location

	#subprocess.call(['cp',target_bed,intermediate_folder])
	#subprocess.call(['cp',param,intermediate_folder])
	subprocess.call(['cp',finalReport_path,sample_folder])
	#subprocess.call(['mv',intermediate_folder+'/tvc_de_novo'+'/processed.bam',sample_folder+'/'+sample+'_'+barcode+'.processed.bam'])
	#subprocess.call(['mv',intermediate_folder+'/tvc_de_novo'+'/processed.bam.bai',sample_folder+'/'+sample+'_'+barcode+'.processed.bam.bai'])
	#subprocess.call(['cp',intermediate_folder+'/_ALAMUT_load_variants.vbs',sample_folder])
	#subprocess.call(['cp',intermediate_folder+'/_ALAMUT_load_processed_bam.vbs',sample_folder])
	#subprocess.call(['cp',intermediate_folder+'/_PRINT_finalReport.vbs',sample_folder])
	
	shutil.make_archive(intermediate_folder,'zip',intermediate_folder)
	#shutil.rmtree(intermediate_folder)

	print " [%s] ... Done" % (time.strftime("%H:%M:%S"))
	
	if options.full_run:
		run_name = run_folder.split('/')[-1]
		if 'HD802-HD748' in sample.upper():
			subprocess.call(['python','/DATA/work/scripts/HD802-HD748_check.py',intermediate_folder+'/format_diag.tsv',sample,run_name])
		elif 'ACROMETRIX' in sample.upper():
			subprocess.call(['python','/DATA/work/scripts/Acrometrix_check.py',finalReport_path,sample,run_name])
		elif run_type == 'SBT':
			subprocess.call(['python','/DATA/work/scripts/collect_variants_MET_intron_13-14.py',finalReport_path,sample,run_name])
		for control_name in control_names:
			if control_name in sample.upper():
				checkconta_bamlist.append(bamfile)

if options.full_run and os.path.isfile(run_folder+'/barcodes.json'):
	
	### VBS scripts for printing
	print " [%s] Making VBS for fast printing..." % (time.strftime("%H:%M:%S"))
	subprocess.call(['python','/DATA/work/finalReport/make_vbs_print.py',run_folder,run_folder+'/barcodes.json'])
	print " [%s] ... Done" % (time.strftime("%H:%M:%S"))
	
	#########################
	## CHECK CONTAMINATION ##
	#########################
	
	checkconta_folder = run_folder+'/_Check-contamination'
	if not os.path.isdir(checkconta_folder):
		subprocess.call(['mkdir', checkconta_folder])
		
	for controlbam in checkconta_bamlist:
		cmd = subprocess.Popen(['python','/DATA/work/checkContamination/Check_contamination.py','--bam', controlbam,'--read-len', '100',], stdout=subprocess.PIPE)
		out, err = cmd.communicate()
		print "\t -Check_contamination <%s> done, ERR: %s" % (controlbam.split('/')[-1].split('_IonXpress')[0],err)
		
	##################
	## PLOTCOVERAGE ##
	##################
	
	cmd = subprocess.Popen(['python', '/DATA/work/plotCoverage/plotCoverage.py','--run-folder', run_folder], stdout=subprocess.PIPE)	
	out, err = cmd.communicate()
	print "\t -plotCoverage done, OUT: %s, ERR: %s"% (out, err)
		
	##########
	## DUPG ##
	##########
	
	cmd = subprocess.Popen(['python', '/DATA/work/dupG/dupG.py','--run-folder', run_folder], stdout=subprocess.PIPE)	
	out, err = cmd.communicate()
	print "\t -dupG done, OUT: %s, ERR: %s"% (out, err)
	
	##############
	## CHECKMUT ##
	##############
	
	subprocess.call(['python','/DATA/work/scripts/checkMut.py','--run-folder',options.full_run,'--chrom','chr17','--start','37880973','--end','37881000','--wt','GAAGCATACGTGATGGCT','--variant','GAAGCATACGTGATGGCATACGTGATGGCT','--mutname','ERBB2_2313_2324dup_COSM20959'])
	subprocess.call(['python','/DATA/work/scripts/checkMut.py','--run-folder',options.full_run,'--chrom','chr3','--start','178952080','--end','178952090','--wt','TGCACATCATG','--variant','TGCACGTCATG','--mutname','PIK3CA_3140_A_G'])
	
	#################
	## VARIANTBASE ##
	#################
	
	cmd = subprocess.Popen(['python', '/DATA/work/variantBase/variantCollector.py','--run-folder', run_folder,'--sbt-wait'], preexec_fn=os.setpgrp)
	print "\t -variantCollector launched"
