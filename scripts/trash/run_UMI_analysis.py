#!/usr/bin/python
import sys
import os
import time
import glob
import json
import subprocess
import zipfile
from optparse import OptionParser

### GATHERING PARAMETERS ############################################################

parser = OptionParser()
parser.add_option('-u', '--ubam',		help="Input bam file for SINGLE BAM ANALYSIS (must contain flow signal from TSS)", dest='ubam')
parser.add_option('-r', '--full-run',	help="Run folder path  for FULL RUN ANALYSIS", dest='full_run') 
parser.add_option('-t', '--run-type',	help="Run type : SBT, LAM, TP53-HEMATO... (OPTIONAL, NOT NEEDED if barcodes.json is present in run folder)", dest='run_type') 
(options, args) = parser.parse_args()

if options.ubam and options.full_run:
	sys.stderr.write("[run_analysis.py] Error: <--bam> and <--full-run> are not compatibles\n")
	sys.exit()
if options.full_run:
	ubamlist = glob.glob(options.full_run+'/*/*.ubam')
	ubamlist = [item for item in ubamlist if not 'processed' in item]
elif options.ubam:
	ubamlist = [options.ubam]
else:
	sys.stderr.write("[run_analysis.py] Error: no <--bam> or <--full-run> specified\n")
	sys.exit()
	
with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)
reference = global_param['ref']
fp_file = global_param['fp_file']

control_names = ['H2O','H20','NTC'] # liste des noms possibles pour les temoins negatifs
checkconta_bamlist = []
sampleask = True

### RUN ANALYSIS ######################################################################

bam_data = {}
	
for ubamfile in sorted(ubamlist) :
	sample = ubamfile.split('/')[-1].split('_IonXpress')[0]
	barcode = 'IonXpress_' + ubamfile.split('IonXpress_')[-1].split('.ubam')[0]
	sample_folder = os.path.dirname(ubamfile)
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
	results_name = ubamfile.split('/')[-1].replace(' ','@') # necessaire, variant_caller_pipeline.py plante si espace dans output_dir
	results_folder = '/DATA/work/results/' + results_name
	if not os.path.isdir(results_folder):
		subprocess.call(['mkdir', results_folder])
		
	## UBAM 2 FASTQ ## 
	
	fastqfile = ubamfile.replace('.ubam','.fastq')
	if not os.path.isfile(fastqfile):
		subprocess.call(['bedtools', 'bamtofastq', '-i', ubamfile, '-fq', fastqfile)
	
	bam_data[ubamfile] = {'barcode':barcode,'sample':sample,'fastq':fastqfile,'run_type':run_type,'sample_folder':sample_folder,'results_folder':results_folder,'target_bed':target_bed,'bed_merged':bed_merged,'param':param,'param_hotspot_only':param_hotspot_only,'hotspot_vcf':hotspot_vcf,'ds_file':ds_file}

		
###########################
## SM_COUNTER PARAM FILE ##
###########################

sm_counter_param_template = open('/DATA/work/reference_files/run_sm_counter_v2.params.iontorrent.txt','r')
sm_counter_param_run = open('/DATA/tests/run_sm_counter_v2.params.iontorrent.txt','r')
for line in sm_counter_param_template:
	sm_counter_param_run.write(line+'\n')

for ubamfile in bam_data:
	sm_counter_param_run.write('[%s]\n' % ubamfile)
	sm_counter_param_run.write('uBam = %s\n' % (bam_data[ubamfile]['fastq']))
	sm_counter_param_run.write('readFile1 = %s\n' % (bam_data[ubamfile]['fastq']))
	sm_counter_param_run.write('readFile2 = \n')
	sm_counter_param_run.write('primerFile = %s\n' % (bam_data[ubamfile]['primerFile']))		# TODO # TODO # TODO
	sm_counter_param_run.write('roiBedFile = %s\n' % (bam_data[ubamfile]['roiBedFile']))		# TODO # TODO # TODO
	sm_counter_param_run.write('platform = IonTorrent\n')
	sm_counter_param_run.write('runCNV = False\n')
	sm_counter_param_run.write('sampleType =  Single\n')
	sm_counter_param_run.write('duplex = False\n')


######################
## LANCEMENT DOCKER ##
######################

container_id = '675af1eab2d1'
for ubamfile in bam_data:
	subprocess.Popen(['docker', 'exec', '-it', container_id, 'python', '/srv/qgen/code/qiaseq-dna/run_qiaseq_dna.py', 'run_sm_counter_v2.params.iontorrent.txt', 'v2', 'single', bam_data[bamfile]['barcode'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)#  > run.log 2>&1
	out, err = cmd.communicate()

	
#######################
## COVERAGE ANALYSIS ## 
#######################

#print " [%s] Coverage Analysis ..." % (time.strftime("%H:%M:%S"))	
#for bamfile in bam_data:
	#coverage_folder = '%s/coverage' % bam_data[bamfile]['results_folder']
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
	##cmd = subprocess.Popen(['python','/DATA/work/CNV/ioncopy/make_ioncopy_input.py','--run-folder', run_folder])
	##out, err = cmd.communicate()
	##cmd = subprocess.call(['bash','/DATA/work/CNV/run_cna.sh', run_folder])
	#cmd = subprocess.call(['python','/DATA/work/CNV/run_cna.py', run_folder])
	
	
#### RESTE DE L'ANALYSE PATIENT PAR PATIENT ###	
#for bamfile in sorted(ubamlist) :
	#print " [%s] Processing %s ..." % (time.strftime("%H:%M:%S"),bamfile.split('/')[-1])
	
	####################################################################
	## FOR SKIPING SAMPLES # WARNING check-contamination will not work #
	##if sampleask:
		##proceed = raw_input('continue? (y/n/stopask)\n')              
		##if proceed == 'y' or proceed == 'yes':                        
			##pass
		##elif proceed == 'stopask':
			##sampleask = False
			##pass                                                       
		##else:                                                         
			##continue                                                  
	####################################################################
	
	#print "\t -run_type =  %s" % run_type	
	#barcode = bam_data[bamfile]['barcode']
	#sample = bam_data[bamfile]['sample']
	#run_type = bam_data[bamfile]['run_type']
	#sample_folder = bam_data[bamfile]['sample_folder']
	#results_folder = bam_data[bamfile]['results_folder']
	#target_bed = bam_data[bamfile]['target_bed']
	#bed_merged = bam_data[bamfile]['bed_merged']
	#param = bam_data[bamfile]['param']
	#param_hotspot_only = bam_data[bamfile]['param_hotspot_only']
	#hotspot_vcf = bam_data[bamfile]['hotspot_vcf']
	#ds_file = bam_data[bamfile]['ds_file']
	
	####################
	### VARIANTCALLER ##
	####################
	
	#if not os.path.isdir(results_folder+'/tvc_de_novo'):
		#subprocess.call(['mkdir', results_folder+'/tvc_de_novo'])
	#if not os.path.isdir(results_folder+'/tvc_only_hotspot'):
		#subprocess.call(['mkdir', results_folder+'/tvc_only_hotspot'])
	
	## RUN 1 : DE NOVO
	
	#cmd = subprocess.Popen([
		#'python', '/DATA/work/variantCaller/bin/variant_caller_pipeline.py',
		#'--input-bam',       bamfile,
		#'--output-dir',      results_folder+'/tvc_de_novo',
		#'--reference-fasta', reference,
		#'--region-bed',      bed_merged,
		#'--primer-trim-bed', target_bed,
		#'--parameters-file', param,
		#'--error-motifs',	'/DATA/work/variantCaller/share/TVC/sse/motifset.txt',
		#'--postprocessed-bam',results_folder+'/tvc_de_novo'+'/processed.bam',
		#], stdout=open(results_folder+'/tvc_de_novo'+'/variant_caller_pipeline.stdout.txt','w'),stderr=open(results_folder+'/tvc_de_novo'+'/variant_caller_pipeline.stderr.txt','w'))
	#out, err = cmd.communicate()
	#print "\t -variant_caller_pipeline.py <de novo> done, OUT: %s, ERR: %s"% (out, err)

	#subprocess.call(['gzip','-d','-c',results_folder+'/tvc_de_novo'+'/TSVC_variants.vcf.gz'],stdout=open(results_folder+'/tvc_de_novo'+'/TSVC_variants.vcf','w'))

	#cmd = subprocess.Popen([
		#'python', '/DATA/work/variantCaller/bin/generate_variant_tables.py',
		#'--input-vcf',		results_folder+'/tvc_de_novo'+'/TSVC_variants.vcf',
		#'--region-bed',		target_bed,
		#'--output-xls',		results_folder+'/tvc_de_novo'+'/output.xls',
		#'--alleles2-xls',	results_folder+'/tvc_de_novo'+'/alleles.xls'
		#], stdout=subprocess.PIPE)
	#out, err = cmd.communicate()
	#print "\t -generate_variant_tables.py <de novo> done, OUT: %s, ERR: %s"%(out, err)	
	
	## RUN 2 : ONLY HOTSPOT
	
	#if hotspot_vcf != '':
		#cmd = subprocess.Popen([
			#'python', '/DATA/work/variantCaller/bin/variant_caller_pipeline.py',
			#'--input-bam',       bamfile,
			#'--output-dir',      results_folder+'/tvc_only_hotspot',
			#'--reference-fasta', reference,
			#'--region-bed',      bed_merged,
			#'--primer-trim-bed', target_bed,
			#'--parameters-file', param_hotspot_only,
			#'--error-motifs',	'/DATA/work/variantCaller/share/TVC/sse/motifset.txt',
			#'--hotspot-vcf',     hotspot_vcf,
			#], stdout=open(results_folder+'/tvc_only_hotspot'+'/variant_caller_pipeline.stdout.txt','w'),stderr=open(results_folder+'/tvc_only_hotspot'+'/variant_caller_pipeline.stderr.txt','w'))
		#out, err = cmd.communicate()
		#print "\t -variant_caller_pipeline.py <only hotspot> done, OUT: %s, ERR: %s"% (out, err)

		#subprocess.call(['gzip','-d','-c',results_folder+'/tvc_only_hotspot'+'/TSVC_variants.vcf.gz'],stdout=open(results_folder+'/tvc_only_hotspot'+'/TSVC_variants.vcf','w'))

		#cmd = subprocess.Popen([
			#'python', '/DATA/work/variantCaller/bin/generate_variant_tables.py',
			#'--input-vcf',		results_folder+'/tvc_only_hotspot'+'/TSVC_variants.vcf',
			#'--region-bed',		target_bed,
			#'--hotspots',
			#'--output-xls',		results_folder+'/tvc_only_hotspot'+'/output.xls',
			#'--alleles2-xls',	results_folder+'/tvc_only_hotspot'+'/alleles.xls'
			#], stdout=subprocess.PIPE)
		#out, err = cmd.communicate()
		#print "\t -generate_variant_tables.py <only hotspot> done, OUT: %s, ERR: %s"%(out, err)	
			
	########################
	### VARIANTANNOTATION ##
	########################

	#cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/generate_annovar_avinput.py',results_folder+'/tvc_de_novo/alleles.xls'], stdout=open(results_folder+'/tvc_de_novo/generate_annovar_avinput.stdout.txt','w'), stderr=open(results_folder+'/tvc_de_novo/generate_annovar_avinput.stderr.txt','w'))
	#out, err = cmd.communicate()
	#print "\t -generate_annovar_avinput.py <de novo> done, OUT: %s, ERR: %s"%(out, err)
	#if hotspot_vcf != '':
		#cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/generate_annovar_avinput.py',results_folder+'/tvc_only_hotspot/alleles.xls'], stdout=open(results_folder+'/tvc_only_hotspot/generate_annovar_avinput.stdout.txt','w'), stderr=open(results_folder+'/tvc_only_hotspot/generate_annovar_avinput.stderr.txt','w'))
		#out, err = cmd.communicate()
		#print "\t -generate_annovar_avinput.py <only hotspot> done, OUT: %s, ERR: %s"%(out, err)
	#cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/merge_annovar_avinput.py',results_folder], stdout=open(results_folder+'/merge_annovar_avinput.stdout.txt','w'), stderr=open(results_folder+'/merge_annovar_avinput.stderr.txt','w'))
	#out, err = cmd.communicate()
	#print "\t -merge_annovar_avinput.py done, OUT: %s, ERR: %s"%(out, err)

	#cmd = subprocess.Popen([
		#'perl', '/DATA/work/variantAnnotation/annovar/table_annovar.pl', results_folder+'/variants.merged.avinput', '/DATA/work/variantAnnotation/annovar/humandb/',
		#'-buildver', 'hg19',
		#'-out', results_folder+'/variants.annotated',
		#'-protocol', 'refGene,cosmic82,avsnp147,intervar_20170202,clinvar_20160302,nci60,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_sas,dbnsfp33a',
		#'-operation', 'g,f,f,f,f,f,f,f,f,f,f,f,f,f',
		#'-argument', '--hgvs,,,,,,,,,,,,,',
		#'-nastring', '.',
		#'-polish',
		#'-xref', '/DATA/work/variantAnnotation/annovar/example/gene_xref.txt',
		#], stdout=open(results_folder+'/table_annovar.stdout.txt','w'), stderr=open(results_folder+'/table_annovar.stderr.txt','w'))
	#out, err = cmd.communicate()
	#print "\t -table_annovar.pl done, ERR: %s"% err

	### format diag
	#cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/format_diag.py',results_folder+'/variants.annotated.hg19_multianno.txt',target_bed], stdout=subprocess.PIPE)
	#out, err = cmd.communicate()
	#print "\t -format_diag.py done, OUT: %s, ERR: %s"%(out, err)
	
	#### OPTIONAL : format VEP with all possible transcripts
	##if run_type in ['Lymphome_B','Lymphome_T','LAM','SBT']:
	#cmd = subprocess.Popen(['python', '/DATA/work/variantAnnotation/VEP/vep.py', results_folder+'/variants.merged.avinput', results_folder+'/format_vep.tsv'], stdout=open(results_folder+'/vep.stdout.txt','w'), stderr=open(results_folder+'/vep.stderr.txt','w'))
	#out, err = cmd.communicate()
	#print "\t -vep.py done, ERR: %s"% err

	##################
	### FINALREPORT ##
	##################
	#if run_type == 'old_SBT_v5':
		#run_type = 'SBT'
	
	#finalreport_name = sample+'_'+barcode+'.finalReport.xlsx'
	#finalReport_path = results_folder+'/'+finalreport_name
	#cmd = subprocess.Popen([
		#'python','/DATA/work/finalReport/finalReport.py',
		#'--input-tsv', results_folder+'/format_diag.tsv',
		#'--output', finalReport_path,
		#'--sample',	sample,
		#'--barcode', barcode,
		#'--sample_folder', sample_folder,
		#'--run_type', run_type,
		#'--target-bed', target_bed,
		#'--fp-file', fp_file,
		#'--ds-file', ds_file,
	#], stdout=open(results_folder+'/finalReport.stdout.txt','w'))
	#out, err = cmd.communicate()
	#print "\t -finalReport.py done, ERR: %s"% err

	#################################################
	## copy finalreport and vbs to bamfile location

	#subprocess.call(['cp',finalReport_path,sample_folder])
	#subprocess.call(['cp',results_folder+'/tvc_de_novo'+'/processed.bam',sample_folder+'/'+sample+'_'+barcode+'.processed.bam'])
	#subprocess.call(['cp',results_folder+'/tvc_de_novo'+'/processed.bam.bai',sample_folder+'/'+sample+'_'+barcode+'.processed.bam.bai'])
	#subprocess.call(['cp',results_folder+'/_ALAMUT_load_variants.vbs',sample_folder])
	#subprocess.call(['cp',results_folder+'/_ALAMUT_load_processed_bam.vbs',sample_folder])
	#subprocess.call(['cp',results_folder+'/_PRINT_finalReport.vbs',sample_folder])
	
	#with zipfile.ZipFile('%s/intermediate_files.zip' % sample_folder,'w') as z:
		#z.write(target_bed,'target.bed')
		#z.write(param,'vc_parameters.json')
		#z.write('%s/coverage/%s_%s.amplicon.cov.xls' % (results_folder,sample,barcode),'amplicon.cov.xls')
		#z.write('%s/format_diag.tsv' % results_folder,'format_diag.tsv')
		#z.write('%s/variants.merged.avinput' % results_folder,'variants.merged.avinput')
		#z.write('%s/tvc_de_novo/alleles.xls' % results_folder,'de_novo_alleles.xls')
		#z.write('%s/tvc_de_novo/TSVC_variants.vcf' % results_folder,'de_novo_TSVC_variants.vcf')
		#if hotspot_vcf != '':
			#z.write(hotspot_vcf,'hotspot.vcf')
			#z.write('%s/tvc_only_hotspot/alleles.xls' % results_folder,'only_hotspot_alleles.xls')
			#z.write('%s/tvc_only_hotspot/TSVC_variants.vcf' % results_folder,'only_hotspot_TSVC_variants.vcf')

	#print " [%s] ... Done" % (time.strftime("%H:%M:%S"))
	
	#if options.full_run:
		#run_name = run_folder.split('/')[-1]
		#if 'HD802-HD748' in sample.upper():
			#subprocess.call(['python','/DATA/work/scripts/HD802-HD748_check.py',results_folder+'/format_diag.tsv',sample,run_name])
		#elif 'ACROMETRIX' in sample.upper():
			#subprocess.call(['python','/DATA/work/scripts/Acrometrix_check.py',finalReport_path,sample,run_name])
		#elif run_type == 'SBT':
			#subprocess.call(['python','/DATA/work/scripts/collect_variants_MET_intron_13-14.py',finalReport_path,sample,run_name])
			#subprocess.call(['python','/DATA/work/scripts/COSM20959.py','--run-folder',options.full_run])
		#for control_name in control_names:
			#if control_name in sample.upper():
				#checkconta_bamlist.append(bamfile)

#if options.full_run and os.path.isfile(run_folder+'/barcodes.json'):
	#print " [%s] Making VBS for fast printing..." % (time.strftime("%H:%M:%S"))
	#subprocess.call(['python','/DATA/work/finalReport/make_vbs_print.py',run_folder,run_folder+'/barcodes.json'])
	#print " [%s] ... Done" % (time.strftime("%H:%M:%S"))
	
	##########################
	### CHECK CONTAMINATION ##
	##########################
	
	#checkconta_folder = run_folder+'/_Check-contamination'
	#if not os.path.isdir(checkconta_folder):
		#subprocess.call(['mkdir', checkconta_folder])
		
	#for controlbam in checkconta_bamlist:
		#cmd = subprocess.Popen(['python','/DATA/work/checkContamination/Check_contamination.py','--bam', controlbam,'--read-len', '100',], stdout=subprocess.PIPE)
		#out, err = cmd.communicate()
		#print "\t -Check_contamination <%s> done, ERR: %s" % (controlbam.split('/')[-1].split('_IonXpress')[0],err)
		
	###################
	### PLOTCOVERAGE ##
	###################
	
	#cmd = subprocess.Popen(['python', '/DATA/work/plotCoverage/plotCoverage.py','--run-folder', run_folder], stdout=subprocess.PIPE)	
	#out, err = cmd.communicate()
	#print "\t -plotCoverage done, OUT: %s, ERR: %s"% (out, err)
		
	###########
	### DUPG ##
	###########
	
	#cmd = subprocess.Popen(['python', '/DATA/work/dupG/dupG.py','--run-folder', run_folder], stdout=subprocess.PIPE)	
	#out, err = cmd.communicate()
	#print "\t -dupG done, OUT: %s, ERR: %s"% (out, err)
	
	##################
	### VARIANTBASE ##
	##################
	
	#cmd = subprocess.Popen(['python', '/DATA/work/variantBase/variantCollector.py','--run-folder', run_folder,'--sbt-wait'], preexec_fn=os.setpgrp)
	#print "\t -variantCollector launched"
