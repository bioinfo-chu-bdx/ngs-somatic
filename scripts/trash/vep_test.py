##!/usr/bin/env python

import requests, sys
 
server = "http://grch37.rest.ensembl.org"
ext = "/vep/human/hgvs"
headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
r = requests.post(server+ext, headers=headers, data='{ "hgvs_notations" : ["7:g.116340262A>G"] }')
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
print repr(decoded)

#import json
#import requests
#import sys
#import csv

#def convert_variation_description(chrm,pos,ref,alt):
	#chrm = chrm.replace('chr','')
	#pos = int(pos)
	#if ref == '-': # insertion       exemple : g.4426_4427insAG
		#variation = "%s:g.%s_%sins%s" % (chrm,str(pos),str(pos+1),alt)
	#elif alt == '-': # deletion       exemple : chr11:g.111959693_111959696delTGAC
		#if len(ref) > 1:
			#variation = "%s:g.%s_%sdel%s" % (chrm,pos,str(pos+(len(ref)-1)),ref)
		#else:
			#variation = "%s:g.%sdel%s" % (chrm,pos,ref)
	#elif len(ref) > 1 or len(alt) > 1: # delins     chr11:g.111959693delinsTGAC  ou chr11:g.111959693_111959694delinsTG
		#if len(ref) > 1:
			#variation = "%s:g.%s_%sdelins%s" % (chrm,pos,str(pos+(len(ref)-1)),alt)
		#else:
			#variation = "%s:g.%sdelins%s" % (chrm,str(pos),alt)
	#else: # SNP
		#variation = "%s:g.%s%s>%s" % (chrm,str(pos),ref,alt)
	#return variation


##vepdata = subprocess.call("wget -q --header='Content-type:application/json' 'http://grch37.rest.ensembl.org/vep/human/hgvs/3:g.38182641T>C?refseq=1'  -O -",shell=True)
##d = json.loads(vepdata)

##for desc in d[0]['transcript_consequences']:
	##refseq_id = desc['transcript_id']
	##if 'NM_002468' in refseq_id:
		##consequence = desc['consequence_terms'][0] # => "missense_variant"

### https://www.ensembl.org/info/genome/variation/predicted_data.html#consequences

#avinput_path = sys.argv[1]
#avinput_file = open(avinput_path,'r')
#avinput_reader = csv.reader(avinput_file,delimiter='\t')

#var_list = []
#for avinput_line in avinput_reader:
	##comm = avinput_line[5]
	##var_freq = comm.split('var_freq=')[-1].split(',')[0]
	##var_cov = comm.split('var_cov=')[-1].split(',')[0]
	##pos_cov = comm.split('pos_cov=')[-1].split(',')[0]
	##var_class = comm.split('var_class=')[-1].split(',')[0]
	##vc_gene = comm.split('gene=')[-1].split(',')[0].split('_ex')[0]
	##vc_call = comm.split('call=')[-1].split(',')[0]
	#chrm = avinput_line[0]
	#pos = avinput_line[1]
	#ref = avinput_line[3]
	#alt = avinput_line[4]

	#var_string = convert_variation_description(chrm,pos,ref,alt)
	#var_list.append(var_string)

#server = "http://rest.ensembl.org"
#ext = "/vep/human/hgvs"
#headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
#r = requests.post(server+ext, headers=headers, data='{ "hgvs_notations" : %s }')
 
#if not r.ok:
  #r.raise_for_status()
  #sys.exit()
 
#decoded = r.json()
#print repr(decoded)
