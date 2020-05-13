#!/usr/bin/env python
import csv
import json
import re
import urllib2

with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)
with open(global_param['NM_data'], 'r') as nmdata:
	NM_version = json.load(nmdata)
sbt_sensitivity = open('/media/stuff/database_hotspot_sensitivity_27-09-2018.updated.tsv','r')
sbt_sensitivity_updated = open('/media/stuff/database_hotspot_sensitivity_01-10-2018.updated.tsv','w')
sbt_sensi_reader = csv.reader(sbt_sensitivity,delimiter='\t')
sbt_sensi_writer = csv.writer(sbt_sensitivity_updated,delimiter='\t')

cosmicMutantExport = open('/media/stuff/COSMIC/CosmicMutantExport.tsv','r')
cosmicMutantExport_reader = csv.reader(cosmicMutantExport,delimiter='\t')

cosmic_lines = []
cosmicMutantExport_reader.next()
for cline in cosmicMutantExport_reader:
	cosmic_lines.append([cline[0],cline[17],cline[16],cline[23]])

header = sbt_sensi_reader.next()
sbt_sensi_writer.writerow(header)

customindex = 1

for line in sbt_sensi_reader:
	nm = line[4]
	gene = line[5]
	c_pos = line[7]
	if c_pos == '':
		print "%s:%s\t%s-%s\t\t%s" % (nm,c_pos,'','','custom'+str(customindex))
		sbt_sensi_writer.writerow([line[0],line[1],line[2],'custom'+str(customindex),line[3],line[4],line[5],line[6],line[7],line[8],line[9]])
		customindex = customindex+1
		continue
	c_pos_pattern = r"%s" % re.escape(c_pos) # snv ou ins ok
	if 'delins' in c_pos:
		before_delins = c_pos.split('delins')[0]
		after_delins = c_pos.split('delins')[-1]
		c_pos_pattern = r"%s[ATGC]+>%s" % (re.escape(before_delins),re.escape(after_delins))
	elif 'del' in c_pos:
		before_del = c_pos.split('del')[0]
		c_pos_pattern = r"%sdel[ATGC]+|%sdel\d+" % (re.escape(before_del),re.escape(before_del))
	
	cosmic_found = []
	start = ''
	stop = ''
	cosmicMutantExport.seek(0)
	cosmicMutantExport_reader.next()
	for cosmic_line in cosmic_lines:
		gene_name = cosmic_line[0]
		#accession_name = cosmic_line[1].split('.')[0]
		Mutation_CDS = cosmic_line[1]	
		mutation_id = cosmic_line[2]
		mutation_genomic_position = cosmic_line[3]

		if not gene_name == gene:
			continue
		if re.match(c_pos_pattern,Mutation_CDS):
			cosmic_found.append(mutation_id)
			start = cosmic_line[3].split(':')[-1].split('-')[0]
			stop = cosmic_line[3].split('-')[-1]		
	cosmic_found = list(set(cosmic_found))
	if cosmic_found :
		print "%s:%s\t%s-%s\t\t%s" % (nm,c_pos,start,stop,cosmic_found)
		sbt_sensi_writer.writerow([line[0],start,stop,','.join(cosmic_found),nm,gene,line[6],c_pos,line[8],line[9]])
	else:
		try:
			for v in range(int(NM_version[nm]['version']),0,-1):
				try:
					nmv = nm + '.' + str(int(v))
					url_mutalyzer = "https://mutalyzer.nl/json/numberConversion?build=hg19;variant=%s:%s" % (nmv,c_pos)
					f = urllib2.urlopen(url_mutalyzer)
					d = json.loads(f.read())
				except:
					continue
				break	
			start = d[0].split('g.')[-1].split('_')[0].split('A')[0].split('T')[0].split('G')[0].split('C')[0].split('del')[0].split('dup')[0].split('inv')[0].split('ins')[0]
			stop = d[0].split('g.')[-1].split('_')[-1].split('A')[0].split('T')[0].split('G')[0].split('C')[0].split('del')[0].split('dup')[0].split('inv')[0].split('ins')[0]
			
			print "%s:%s\t%s-%s\t\t%s" % (nm,c_pos,start,stop,'custom'+str(customindex))
			sbt_sensi_writer.writerow([line[0],start,stop,'custom'+str(customindex),nm,gene,line[6],c_pos,line[8],line[9]])
			customindex = customindex+1
		except:
			print "%s:%s\t%s-%s\t\t%s" % (nm,c_pos,'','','custom'+str(customindex))
			sbt_sensi_writer.writerow([line[0],'','','custom'+str(customindex),nm,gene,line[6],c_pos,line[8],line[9]])
			customindex = customindex+1
		
