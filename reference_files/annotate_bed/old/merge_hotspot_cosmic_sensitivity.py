#!/usr/bin/env python
import csv
import json
import re
import urllib2
import subprocess

sbt_sensitivity = open('/media/stuff/database_hotspot_sensitivity_01-10-2018.updated.tsv','r')
sbt_sensi_reader = csv.reader(sbt_sensitivity,delimiter='\t')

cosmic_hotspot_generated = open('/DATA/work/reference_files/Annotate_bed_script/SBT_results3_hotspots.vcf','r')
cosmic_hotspot_generated_reader = csv.reader(cosmic_hotspot_generated,delimiter='\t')

merged_final_hotspot = open('/DATA/work/reference_files/Annotate_bed_script/SBT_final_hotspots_merged2.vcf','w')
merged_final_hotspot_writer = csv.writer(merged_final_hotspot,delimiter = '\t')

lineindex = 0

sensi2write = []

sbt_sensi_reader.next()
for line in sbt_sensi_reader:
	lineindex = lineindex + 1
	chrom = line[0].replace('chr','')
	pos = line[1]
	mut_id = line[3]
	gene = line[5]
	sensi = line[9]
	if sensi in ['exclu','non incluable','exclu protocole ACSE','polymorphisme']:
		continue
	strand = '+'
	reverse = False
	if gene in ['KRAS','NRAS','BRAF','HRAS','POLE','ALK','FGFR2']:
		strand = '-'
		reverse = True
	cpos = line[7]
	ppos = line[8]
	ref = line[2].split('/')[0]
	alt = line[2].split('/')[-1]
	if (ref == '') or (alt == '') or (pos == ''):
		print "- line %s dismissed (not enough data)" % lineindex
		continue
	if reverse :
		ref = ref.replace('A','W').replace('T','X').replace('G','Y').replace('C','Z')
		alt = alt.replace('A','W').replace('T','X').replace('G','Y').replace('C','Z')
		ref = ref.replace('W','T').replace('X','A').replace('Y','C').replace('Z','G')
		alt = alt.replace('W','T').replace('X','A').replace('Y','C').replace('Z','G')
		ref = ref[::-1]
		alt = alt[::-1]
	# FIND ANCHOR if '-'
	if ref == '-':
		faidx = subprocess.check_output(['samtools','faidx','/DATA/work/hg19.fasta','%s:%s-%s' % (line[0],int(pos)-1,int(pos)-1)])
		anchor = faidx.split('\n')[-2]
		ref = anchor
		alt = anchor + alt
		pos = str(int(pos)-1)
	if alt == '-':
		faidx = subprocess.check_output(['samtools','faidx','/DATA/work/hg19.fasta','%s:%s-%s' % (line[0],int(pos)-1,int(pos)-1)])
		anchor = faidx.split('\n')[-2]
		alt = anchor
		ref = anchor + ref
		pos = str(int(pos)-1)
	
	sensi2write.append([chrom,pos,mut_id,ref,alt,'.','.','GENE=%s;STRAND=%s,CDS=%s;AA=%s' % (gene,strand,cpos,ppos)])
		
for line in cosmic_hotspot_generated_reader:
	this = (line[0],line[1],line[3],line[4])
	found = False
	for sensi_line in sensi2write:
		that = (sensi_line[0],sensi_line[1],sensi_line[3],sensi_line[4])
		if this == that:
			#print "found"
			found = True
			continue
	if not found :
		print "%s new hotspot, not in sensi file" % str(this)
		sensi2write.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7]])

for line in sensi2write :
	merged_final_hotspot_writer.writerow(line)
