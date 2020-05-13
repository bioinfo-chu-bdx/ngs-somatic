#!/usr/bin/env python
import subprocess
import sys
import csv
import os

def dot_if_empty(data):
	if data == '-':
		data = '.'
	return data

vep_input = sys.argv[1]
vep_input_file = open(vep_input,'r')
vep_input_reader = csv.reader(vep_input_file,delimiter='\t')

vep_output = os.path.dirname(vep_input)+'/vep_output.tsv'
format_vep = open(os.path.dirname(vep_input)+'/format_vep.tsv','w')

header = ['Comm.Bio.','Chr','Transcript','Gene','Exon','Intron','c.','p.','Consequence','Impact','Var.Freq.','Var.Cov.','Pos.Cov.','Class.','Start.Pos',
'Ref.Seq','Var.Seq','COSMIC','dbSNP','ClinVar.significance','SIFT','PolyPhen','AFR_AF','AMR_AF','EAS_AF','EUR_AF','SAS_AF','AA_AF','EA_AF','VC_CALL']
format_vep.write('\t'.join(header)+'\n')


print "\t -launching vep (offline mode)"
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
	'--variant_class',
	'--sift', 'p',
	'--polyphen', 'p',
	'--af_1kg',
	#'--af_gnomad',
	#'--af_exac',
	'--freq_pop', '1KG_ALL',
	'--freq_pop', 'gnomAD',
	'--af_esp',
	'--symbol',
	'--tab',
	#'--canonical',
	'--pubmed',
	'--fasta', '/DATA/work/variantAnnotation/VEP/cache/homo_sapiens/96_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz',
	'--exclude_predicted', # RETIRE LES TRANSCRITS PREDITS "XM" ou "XR"
	'--no_stats',
	'--verbose',
	
	'-i', vep_input,
	'-o', vep_output

	], stdout=subprocess.PIPE)
#	'--fields', 'Uploaded_variation,Feature,SYMBOL,EXON,INTRON,HGVSc,HGVSp,Consequence,IMPACT',
#	'--format', 'ensembl',
out, err = cmd.communicate()
print "\t -vep done, OUT: %s, ERR: %s"%(out, err)

## FORMAT VEP (FORMAT DIAG VERSION VEP)

# Gathering comments from input file, as they are not keeped in vep output
var_data = {}
for input_line in vep_input_reader:
	comm = input_line[5]
	var_freq = comm.split('var_freq=')[-1].split(',')[0]
	var_cov = comm.split('var_cov=')[-1].split(',')[0]
	pos_cov = comm.split('pos_cov=')[-1].split(',')[0]
	var_class = comm.split('var_class=')[-1].split(',')[0]
	vc_gene = comm.split('gene=')[-1].split(',')[0].split('_ex')[0]
	vc_call = comm.split('call=')[-1].split(',')[0]
	variant_description = ('chr'+input_line[0],input_line[1],input_line[2],input_line[3].split('/')[0],input_line[3].split('/')[1]) # (chr,start,end,ref,alt)
	var_data[variant_description] = {'var_freq':var_freq,'var_cov':var_cov,'pos_cov':pos_cov,'var_class':var_class,'vc_gene':vc_gene,'vc_call':vc_call}

vep_output_file = open(vep_output,'r')
vep_output_reader = csv.reader(vep_output_file,delimiter='\t')

index = {}
for line in vep_output_reader:
	if line[0].startswith('##'):
		continue
	if line[0].startswith('#'): # header
		for i in range(len(line)):
			index[line[i]] = i 
		continue

	ref_nm = line[index['Feature']]
	if not ref_nm.startswith('NM_'):
		continue
			
	comm = '.'
	chrom = 'chr' + line[index['Location']].split(':')[0]
	var_class = line[index['VARIANT_CLASS']]
	if '-' in line[index['Location']]:
		if var_class == 'insertion':
			pos = line[index['Location']].split(':')[1].split('-')[1]
			end = line[index['Location']].split(':')[1].split('-')[0]
		else:
			pos = line[index['Location']].split(':')[1].split('-')[0]
			end = line[index['Location']].split(':')[1].split('-')[1]
	else:
		pos = line[index['Location']].split(':')[1]
		end = pos
	ref = line[index['GIVEN_REF']]
	alt = line[index['Allele']]
	if ref != line[index['USED_REF']]:
		comm = 'Reference error'
	
	#DATA from comments in vep input
	var_freq = var_data[(chrom,pos,end,ref,alt)]['var_freq'] #int(round(float(var_data[(chrom,pos,end,ref,alt)]['var_freq'])))
	var_cov = var_data[(chrom,pos,end,ref,alt)]['var_cov']
	pos_cov = var_data[(chrom,pos,end,ref,alt)]['pos_cov']
	#var_class = line[index['#Uploaded_variation']].split('var_class=')[-1].split(',')[0]
	#vc_gene = var_data[(chrom,pos,end,ref,alt)]['vc_gene']
	vc_call = var_data[(chrom,pos,end,ref,alt)]['vc_call']
	
	#PARSE JSON DATA
	existing_variation = line[index['Existing_variation']].split(',')
	dbSNP = '.'
	COSMIC = '.'
	for item in existing_variation:
		if item.startswith('rs'):
			dbSNP = item
		elif item.startswith('COSM'):
			COSMIC = item
	ClinSig = dot_if_empty(line[index['CLIN_SIG']])
	AFR_AF = dot_if_empty(line[index['AFR_AF']])
	AMR_AF = dot_if_empty(line[index['AMR_AF']])
	EAS_AF = dot_if_empty(line[index['EAS_AF']])
	EUR_AF = dot_if_empty(line[index['EUR_AF']])
	SAS_AF = dot_if_empty(line[index['SAS_AF']])
	AA_AF = dot_if_empty(line[index['AA_AF']])
	EA_AF = dot_if_empty(line[index['EA_AF']])
	gene = line[index['SYMBOL']]
	exon = dot_if_empty(line[index['EXON']].split('/')[0])
	intron = dot_if_empty(line[index['INTRON']].split('/')[0])
	c = dot_if_empty(line[index['HGVSc']].split(':')[-1])
	p = dot_if_empty(line[index['HGVSp']].split(':')[-1])
	consequence = line[index['Consequence']]
	impact = line[index['IMPACT']]
	sift = dot_if_empty(line[index['SIFT']])
	polyphen = dot_if_empty(line[index['PolyPhen']])
	
	format_vep_line = '\t'.join(str(x) for x in [comm,chrom,ref_nm,gene,exon,intron,c,p,consequence,impact,var_freq,var_cov,pos_cov,var_class,pos,ref,alt,COSMIC,dbSNP,ClinSig,sift,polyphen,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,AA_AF,EA_AF,vc_call])
	format_vep.write(format_vep_line+'\n')
	
vep_input_file.close()
format_vep.close()
vep_output_file.close()
