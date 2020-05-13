#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import json
import csv
from optparse import OptionParser
import hgvs.dataproviders.uta
    
parser = OptionParser()
parser.add_option('-b', '--bed',		help="Add new NM from annotated bed", dest='bed')
parser.add_option('-u', '--update',		help="Update NM_data json : search new versions if existant", action='store_true' , default=False ,dest='update')
(options, args) = parser.parse_args()

if not options.bed and not options.update:
	print "no options specified (--bed or/and --update)"
	exit()

with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)
with open(global_param['NM_data'], 'r') as nmdata:
	NM_data = json.load(nmdata)
with open(global_param['NC_data'], 'r') as ncdata:
	NC_data = json.load(ncdata)
	
hdp = hgvs.dataproviders.uta.connect()
	
if options.bed:
	bed_data = {}
	bed_file = open(options.bed,'r')
	bed_lines = bed_file.readlines()
	for i in range(1,len(bed_lines)):
		bed_line = bed_lines[i].split('\t')
		chromosome = bed_line[0]
		nc = NC_data[chromosome]
		bed_nm = bed_line[7].split('TRANSCRIPT=')[-1].split(';')[0].replace('\n','')
		bed_gene = bed_line[7].split('GENE=')[-1].split(';')[0].replace('\n','')
		bed_gene_strand = bed_line[7].split('STRAND=')[-1].split(';')[0].replace('\n','')
		bed_data[bed_gene] = {'nm':bed_nm,'strand':bed_gene_strand}
	
	for gene in bed_data:
		nm = bed_data[gene]['nm']
		if gene in NM_data:
			print "** Warning : Ignoring %s. Gene %s has already %s as default." % (nm, gene, NM_data[gene]['transcript'])
		else:
			NM_data[gene] = {}
			print "* new Gene and NM : %s:%s" % (gene,nm)
			gene_tx = hdp.get_tx_for_gene(gene)
			v = 0
			for item in gene_tx:
				if (item[3].split('.')[0] == nm) and (item[4] == nc):
					v = max(v,int(item[3].split('.')[-1]))			
			NM_data[gene]['transcript'] = nm
			NM_data[gene]['transcriptVersion'] = int(v)
			NM_data[gene]['strand'] = strand
			NM_data[gene]['NCtranscript'] = nc
				
if options.update : 
	for gene in NM_data:
		nm = NM_data[gene]['transcript']
		nc = NM_data[gene]['NCtranscript']
		print "-%s (%s)" % (gene,nm)
		
		gene_tx = hdp.get_tx_for_gene(gene)
		v = 0
		for item in gene_tx:
			if (item[3].split('.')[0] == nm) and (item[4] == nc):
				v = max(v,int(item[3].split('.')[-1]))	
		#
		#
		#CONTINUER ICI
		#
		#
		#
		if str(NM_data[nm]['version']) != str(version_name_checker):
			print "* new version for " + nm
			print "\t %s -> %s" % (NM_data[nm]['version'],version_name_checker)
			NM_data[nm]['version'] = int(version_name_checker)
		for version_position_converter in range(int(version_name_checker),0,-1):
			try: 
				url_mutalyzer_conversion = "https://mutalyzer.nl/json/numberConversion?build=hg19;variant=%s.%s:c.10G>C" % (nm,version_position_converter)
				f = urllib2.urlopen(url_mutalyzer_conversion)
				d = json.loads(f.read())
				start1 = int(d[0].split('g.')[-1].split('_')[0].split('A')[0].split('T')[0].split('G')[0].split('C')[0].split('del')[0].split('dup')[0])
			except:
				continue
			break
		if str(NM_data[nm]['version_position_converter']) != str(version_position_converter):
			print "* new version for %s (position converter) " % nm
			print "\t %s -> %s" % (NM_data[nm]['version_position_converter'],version_position_converter)
			NM_data[nm]['version_position_converter'] = version_position_converter

jsn = open(global_param['NM_data'],'w')
json_text = json.dumps(NM_data,indent=4,sort_keys=True)
jsn.write(json_text)
jsn.close()
