#!/usr/bin/env python
import requests
import sys
import json
import csv
import os
import time
 
def convert_variation_description(chrm,pos,ref,alt):
	chrm = chrm.replace('chr','')
	pos = int(pos)
	if ref == '-': # insertion       exemple : g.4426_4427insAG
		variation = "%s:g.%s_%sins%s" % (chrm,str(pos),str(pos+1),alt)
	elif alt == '-': # deletion       exemple : chr11:g.111959693_111959696delTGAC
		if len(ref) > 1:
			variation = "%s:g.%s_%sdel%s" % (chrm,pos,str(pos+(len(ref)-1)),ref)
		else:
			variation = "%s:g.%sdel%s" % (chrm,pos,ref)
	elif len(ref) > 1 or len(alt) > 1: # delins     chr11:g.111959693delinsTGAC  ou chr11:g.111959693_111959694delinsTG
		if len(ref) > 1:
			variation = "%s:g.%s_%sdelins%s" % (chrm,pos,str(pos+(len(ref)-1)),alt)
		else:
			variation = "%s:g.%sdelins%s" % (chrm,str(pos),alt)
	else: # SNP
		variation = "%s:g.%s%s>%s" % (chrm,str(pos),ref,alt)
	return variation 
 
server = "http://grch37.rest.ensembl.org"

avinput_path = sys.argv[1]
output_path = sys.argv[2]
avinput_file = open(avinput_path,'r')
avinput_reader = csv.reader(avinput_file,delimiter='\t')

isDone = {}
resfile = open(output_path,'w')
jsonfile = open(os.path.dirname(output_path)+'/vep_data.json','w')
header = ['Comm.Bio.','Chr','Transcript','Gene','Exon','Intron','c.','p.','Consequence','Impact','Var.Freq.','Var.Cov.','Pos.Cov.','Class.','Start.Pos',
'Ref.Seq','Var.Seq','COSMIC','dbSNP','ClinVar.significance','SIFT','PolyPhen','AFR_AF','AMR_AF','EAS_AF','EUR_AF','SAS_AF','AA_AF','EA_AF','VC_CALL']
resfile.write('\t'.join(header)+'\n')

#max_request_per_second = 15 # to avoid rate-limit
#request_count = 0

for avinput_line in avinput_reader:
	chrm = avinput_line[0]
	pos = avinput_line[1]
	ref = avinput_line[3]
	alt = avinput_line[4]
	var_string = convert_variation_description(chrm,pos,ref,alt)
	isDone[var_string] = []
	
	comm = avinput_line[5]
	var_freq = comm.split('var_freq=')[-1].split(',')[0]
	var_cov = comm.split('var_cov=')[-1].split(',')[0]
	pos_cov = comm.split('pos_cov=')[-1].split(',')[0]
	var_class = comm.split('var_class=')[-1].split(',')[0]
	vc_gene = comm.split('gene=')[-1].split(',')[0].split('_ex')[0]
	vc_call = comm.split('call=')[-1].split(',')[0]
	ext = "/vep/human/hgvs/%s?refseq=1&numbers=1&hgvs=1" % var_string
	status_code = 0
	print var_string
	while status_code != 200:
		try:
			r = requests.get(server+ext, headers={ "Content-Type" : "application/json"},timeout=5)
		except:
			continue
		status_code = r.status_code
		print status_code
		print r.headers
		if status_code != 200:
			time.sleep(1)
	data = r.json()
	json_text = json.dumps(data, indent=4)
	jsonfile.write(json_text)
	
	#PARSE JSON DATA
	COSMIC = '.'
	dbSNP = '.'
	ClinSig = '.'
	#AF = '.'
	AFR_AF = '.'
	AMR_AF = '.'
	EAS_AF = '.'
	EUR_AF = '.'
	SAS_AF = '.'
	AA_AF = '.'
	EA_AF = '.'
	if 'colocated_variants' in data[0]:
		for cv in data[0]['colocated_variants']:
			if cv['id'].startswith('COSM'):
				COSMIC = cv['id']
			elif cv['id'].startswith('rs'):
				dbSNP = cv['id']
				if 'clin_sig' in cv:
					ClinSig = cv['clin_sig']
				#if 'maf' in cv:
					#AF = cv['maf']
				if 'afr_maf' in cv:
					AFR_AF = cv['afr_maf']
				if 'amr_maf' in cv:
					AMR_AF = cv['amr_maf']
				if 'eas_maf' in cv:
					EAS_AF = cv['eas_maf']
				if 'eur_maf' in cv:
					EUR_AF = cv['eur_maf']
				if 'sas_maf' in cv:
					SAS_AF = cv['sas_maf']
				if 'aa_maf' in cv:
					AA_AF = cv['aa_maf']
				if 'eu_maf' in cv:
					EA_AF = cv['eu_maf']
	if 'transcript_consequences' in data[0].keys():
		for tc in data[0]['transcript_consequences']:
			CommBio = '.'
			RefNM = tc['transcript_id']
			if RefNM in isDone[var_string]:
				continue
			if not RefNM.startswith('NM_'):
				continue
			gene = tc['gene_symbol']
			if 'exon' in tc:
				exon = tc['exon'].split('/')[0]
			else:
				exon = '.'
			if 'intron' in tc:
				intron = tc['intron'].split('/')[0]
			else:
				intron = '.'
			if 'hgvsc' in tc:
				c = tc['hgvsc'].split(':')[-1]
			else:
				c = '.'
			if 'hgvsp' in tc:
				p = tc['hgvsp'].split(':')[-1]
			else:
				p = '.'
			if 'sift_prediction' in tc:
				sift = tc['sift_prediction']
			else:
				sift = '.'
			if 'polyphen_prediction' in tc:
				polyphen = tc['polyphen_prediction']
			else:
				polyphen = '.'
			consequence = ','.join(tc['consequence_terms'])
			impact = tc['impact']
			blabla = '\t'.join(str(x) for x in [CommBio,chrm,RefNM,gene,exon,intron,c,p,consequence,impact,var_freq,var_cov,pos_cov,var_class,pos,ref,alt,COSMIC,dbSNP,ClinSig,sift,polyphen,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,AA_AF,EA_AF,vc_call])
			resfile.write(blabla+'\n')	
			isDone[var_string].append(RefNM)
	else:
		blabla = '\t'.join(str(x) for x in ['no data (bad ref?)',chrm,'.','.','.','.','.','.','.','.',var_freq,var_cov,pos_cov,var_class,pos,ref,alt,'.','.','.','.','.','.','.','.','.','.','.','.',vc_call])
		resfile.write(blabla+'\n')	
			
	#request_count = request_count + 1
	#if request_count == max_request_per_second:
		#request_count = 0
		#time.sleep(1)
			
resfile.close()
jsonfile.close()
