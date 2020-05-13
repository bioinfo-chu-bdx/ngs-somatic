#!/usr/bin/python
import sys
import os
import time
import glob
import json
import subprocess
from optparse import OptionParser

### GATHERING PARAMETERS ############################################################

parser = OptionParser()
parser.add_option('-b', '--bam',		help="Input bam file for SINGLE BAM ANALYSIS (must contain flow signal from TSS)", dest='bam')
parser.add_option('-r', '--run-folder',	help="Run folder path  for FULL RUN ANALYSIS", dest='run_folder') 
parser.add_option('-t', '--run-type',	help="Run type : SBT, LAM, TP53-HEMATO... (OPTIONAL, NOT NEEDED if barcodes.json is present in run folder)", dest='run_type') 
(options, args) = parser.parse_args()

if options.bam and options.run_folder:
	sys.stderr.write("[run_analysis.py] Error: <--bam> and <--full-run> are not compatibles\n")
	sys.exit()
if options.run_folder:
	bamlist = glob.glob(options.run_folder+'/*/*.bam')
	bamlist = [item for item in bamlist if not 'processed' in item]
elif options.bam:
	bamlist = [options.bam]
else:
	sys.stderr.write("[run_analysis.py] Error: no <--bam> or <--full-run> specified\n")
	sys.exit()

with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)

### RUN ANALYSIS ######################################################################

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

	# dossier intermediates files
	intermediate_folder = sample_folder + '/intermediate_files'
	if not os.path.isdir(intermediate_folder):
		subprocess.call(['mkdir', intermediate_folder])
	
###            __              ___     __                         __  ###
### \  /  /\  |__) |  /\  |\ |  |     /  `  /\  |    |    | |\ | / _` ###
###  \/  /~~\ |  \ | /~~\ | \|  |     \__, /~~\ |___ |___ | | \| \__> ###
                                                                  
	print " [%s] Processing %s ..." % (time.strftime("%H:%M:%S"),bamfile.split('/')[-1])
	print "\t * panel %s" % run_type
	
	if not os.path.isdir(intermediate_folder+'/mutect2'):
		subprocess.call(['mkdir', intermediate_folder+'/mutect2'])

	#print "\t - [%s] mutect2 ... (no chunks)" % (time.strftime("%H:%M:%S"))
	#cmd = subprocess.Popen([
		#'/DATA/work/gatk/gatk', 'Mutect2',
		#'-I',	bamfile,
		#'-R',	reference,
		#'-L',	'/DATA/work/reference_files/mutect/Target_ColonLung_v10.intervals',
		#'--max-mnp-distance', '0',
		#'--max-reads-per-alignment-start', '0',
		#'--min-base-quality-score', '1',
		#'--disable-read-filter', 'NotDuplicateReadFilter',
		#'--germline-resource',	'/DATA/work/reference_files/mutect/af-only-gnomad.raw.sites.b37.vcf.gz',
		#'-O', '%s/mutect2/%s_%s.mutect2.vcf.gz' % (intermediate_folder,sample,barcode)
		##'--panel-of-normals', '/DATA/work/reference_files/mutect/colon_lung_pon.vcf.gz'
		#], stdout=open('%s/mutect2/mutect2.stdout.txt'%intermediate_folder,'w'),stderr=open('%s/mutect2/mutect2.stderr.txt'%intermediate_folder,'w'))
	#cmd.communicate()
	
	print "\t - [%s] mutect2 ... (chunks)" % (time.strftime("%H:%M:%S"))
	cmd_list = []
	vcf_chunk_list = []
	for i in range(1,11):
		#interval_chunk = target_bed.replace('reference_files/','reference_files/mutect/').replace('.bed','.chunk%s.intervals'%i)
		interval_chunk = '/DATA/work/reference_files/mutect/colon_lung_splitintervals/000%s-scattered.interval_list'%(i-1)
		vcf_chunk = '%s/mutect2/%s_%s.mutect2.chunk%s.vcf.gz' % (intermediate_folder,sample,barcode,i)
		vcf_chunk_list.append(vcf_chunk)
		print "\t\t - [%s] chunck%s ..." % (time.strftime("%H:%M:%S"),i)
		cmd = subprocess.Popen([
			'/DATA/work/gatk/gatk', 'Mutect2',
			'-I',	bamfile,
			'-R',	reference,
			'-L',	interval_chunk,
			'--max-mnp-distance', '0',
			'--max-reads-per-alignment-start', '0',
			'--min-base-quality-score', '1',
			'--disable-read-filter', 'NotDuplicateReadFilter',
			'--germline-resource',	'/DATA/work/reference_files/mutect/af-only-gnomad.raw.sites.b37.vcf.gz',
			'-O', vcf_chunk
			#'--panel-of-normals', '/DATA/work/reference_files/mutect/colon_lung_pon.vcf.gz'
			], stdout=open('%s/mutect2/mutect2.chunk%s.stdout.txt'%(intermediate_folder,i),'w'),stderr=open('%s/mutect2/mutect2.chunk%s.stderr.txt'%(intermediate_folder,i),'w'))
		cmd_list.append(cmd)
		
	for cmd in cmd_list:
		cmd.wait()
		
	cmd = ['/DATA/work/gatk/gatk', 'GatherVcfs', '-O', '%s/mutect2/%s_%s.mutect2.vcf.gz' % (intermediate_folder,sample,barcode)]
	for vcf_chunk in vcf_chunk_list:
		cmd.append('-I')
		cmd.append(vcf_chunk)
	print "\t\t - [%s] GatherVcfs ..." % time.strftime("%H:%M:%S")
	subprocess.call(cmd,stdout=open(intermediate_folder+'/mutect2/gather_vcf.stdout.txt','w'),stderr=open(intermediate_folder+'/mutect2/gather_vcf.stderr.txt','w'))
	
	cmd = ['/DATA/work/gatk/gatk', 'MergeMutectStats', '-O', '%s/mutect2/%s_%s.mutect2.vcf.gz.stats' % (intermediate_folder,sample,barcode)]
	for vcf_chunk in vcf_chunk_list:
		cmd.append('-stats')
		cmd.append(vcf_chunk+'.stats')
	print "\t\t - [%s] MergeMutectStats ..." % time.strftime("%H:%M:%S")
	subprocess.call(cmd,stdout=open(intermediate_folder+'/mutect2/merge_mutect_stats.stdout.txt','w'),stderr=open(intermediate_folder+'/mutect2/merge_mutect_stats.stderr.txt','w'))
	
	print "\t\t - [%s] Tabix ..." % time.strftime("%H:%M:%S")
	subprocess.call(['tabix','-p','vcf','%s/mutect2/%s_%s.mutect2.vcf.gz' % (intermediate_folder,sample,barcode)],stdout=open(intermediate_folder+'/mutect2/tabix.stdout.txt','w'),stderr=open(intermediate_folder+'/mutect2/tabix.stderr.txt','w'))
	
	print "\t\t - [%s] FilterMutectCalls ..." % time.strftime("%H:%M:%S")
	cmd = subprocess.Popen([
	'/DATA/work/gatk/gatk', 'FilterMutectCalls',
	'-V',	'%s/mutect2/%s_%s.mutect2.vcf.gz' % (intermediate_folder,sample,barcode),
	'-R',	reference,
	'--max-events-in-region', '10',
	'--min-allele-fraction', '0.01',
	'-O', '%s/mutect2/%s_%s.mutect2.filtered.vcf' % (intermediate_folder,sample,barcode)
	#'--panel-of-normals', '/DATA/work/reference_files/mutect/colon_lung_pon.vcf.gz'
	], stdout=open(intermediate_folder+'/mutect2/filter_mutect_calls.stdout.txt','w'),stderr=open(intermediate_folder+'/mutect2/filter_mutect_calls.stderr.txt','w'))
	cmd.communicate()
	
###                 __  ___      ___    __       ###
###  /\  |\ | |\ | /  \  |   /\   |  | /  \ |\ | ###
### /~~\ | \| | \| \__/  |  /~~\  |  | \__/ | \| ###

#####                 __             __  #####
#####  /\  |\ | |\ | /  \ \  /  /\  |__) #####
##### /~~\ | \| | \| \__/  \/  /~~\ |  \ #####
	if not os.path.isdir(intermediate_folder+'/annovar'):
		subprocess.call(['mkdir', intermediate_folder+'/annovar'])
		
	subprocess.call(["awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == \"PASS\") print}' %s/mutect2/%s_%s.mutect2.filtered.vcf > %s/mutect2/%s_%s.mutect2.filtered.pass.vcf"%(intermediate_folder,sample,barcode,intermediate_folder,sample,barcode)], shell=True)
	subprocess.call(['perl','/DATA/work/variantAnnotation/annovar/convert2annovar.pl','-format','vcf4','%s/mutect2/%s_%s.mutect2.filtered.pass.vcf' % (intermediate_folder,sample,barcode),'-outfile','%s/annovar/annovar_intermediate_input.tsv' % intermediate_folder,'-includeinfo'])
	print "\t - [%s] generate_annovar_input_mutect.py ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/generate_annovar_input_mutect.py',
	'--input', '%s/annovar/annovar_intermediate_input.tsv' % intermediate_folder,
	'--output', '%s/annovar/annovar_input.tsv' % intermediate_folder,
	'--run-type', options.run_type
	], stdout=open(intermediate_folder+'/annovar/generate_annovar_input_mutect.stdout.txt','w'), stderr=open(intermediate_folder+'/annovar/generate_annovar_input_mutect.stderr.txt','w'))
	cmd.communicate()

	print "\t - [%s] table_annovar.pl ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen([
		'perl', '/DATA/work/variantAnnotation/annovar/table_annovar.pl', 
		intermediate_folder+'/annovar/annovar_input.tsv', '/DATA/work/variantAnnotation/annovar/humandb/',
		'-buildver', 'hg19',
		'-out', intermediate_folder+'/annovar/variants.annotated',
		'-protocol', 'refGeneWithVer,cosmic89,avsnp150,intervar_20180118,clinvar_20190305,nci60,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eur,gnomad211_genome,exac03,dbnsfp33a',
		'-operation', 'g,f,f,f,f,f,f,f,f,f,f,f',
		'-argument', '--hgvs,,,,,,,,,,,',
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

	print "\t - [%s] generate_vep_input_mutect.py ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen([
		'python', '/DATA/work/variantAnnotation/generate_vep_input_mutect.py',
		'--input', '%s/annovar/annovar_input.tsv' % intermediate_folder,
		'--output', '%s/vep/vep_input.tsv' % intermediate_folder
		], stdout=open(intermediate_folder+'/vep/generate_vep_input.stdout.txt','w'), stderr=open(intermediate_folder+'/vep/generate_vep_input.stderr.txt','w'))
	cmd.communicate()

	print "\t - [%s] run_vep ..." % (time.strftime("%H:%M:%S"))
	#cmd = subprocess.Popen(['python', '/DATA/work/variantAnnotation/VEP/run_vep.py', intermediate_folder+'/vep/vep_input.tsv'], stdout=open(intermediate_folder+'/vep/vep.stdout.txt','w'), stderr=open(intermediate_folder+'/vep/vep.stderr.txt','w'))
	cmd = subprocess.Popen([
		'perl', '/DATA/work/variantAnnotation/VEP/ensembl-vep/vep',
		'--format', 'ensembl',
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
	cmd = subprocess.Popen(['python','/DATA/work/variantAnnotation/format_diag.py',
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
	finalreport_name = sample+'_'+barcode+'.finalReport.xlsx'
	finalReport_path = intermediate_folder+'/finalReport/'+finalreport_name
	
	xmin = '300' # for coverage analysis sheet in finalreport
		
	print "\t - [%s] finalReport.py ..." % (time.strftime("%H:%M:%S"))
	cmd = subprocess.Popen([
		'python','/DATA/work/finalReport/finalReport.py',
		'--bam', bamfile,
		'--run-type', run_type,
		'--output', finalReport_path,
		'--xmin', xmin,
	], stdout=open(intermediate_folder+'/finalReport/finalReport.stdout.txt','w'), stderr=open(intermediate_folder+'/finalReport/finalReport.stderr.txt','w'))
	cmd.communicate()
