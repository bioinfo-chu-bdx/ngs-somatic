#!/usr/bin/python
import os
import sys
import csv
import json
import sqlite3
import subprocess
import uuid
import hgvs.dataproviders.uta

# THIS SCRIPT IS DESIGNED TO FILL AN EMPTY VARIANT BASE WITH PANEL TABLE AND TARGETEDREGION TABLE
# UPDATE LIST UNDER FOR NEW PANELS

# USAGE : python insert_all_panel_variantbase.py

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))
	
with open(global_param['NC_data'], 'r') as ncdata:
	NC_data = json.load(ncdata)
	
refGene_file = open(global_param['RefSeq'],'r')
refGene_reader = csv.reader(refGene_file,delimiter='\t')
hdp = hgvs.dataproviders.uta.connect()
db_path = global_param['VariantBase']

# LIST OF PANEL SORTED BY NEWER > OLDER. (order is important so that old annotations are not written first)
panels = [
	('%s/reference_files/reannoted/Target_ColonLung_v10_IAD172906_231.bed' % pipeline_folder		,'SBT',			10),
	('%s/reference_files/reannoted/LAM2018_IAD143291_236_Designed_with_NM.bed' % pipeline_folder	,'LAM',			4),
	('%s/reference_files/reannoted/TP53.20140108.designed_with_NM.bed' % pipeline_folder			,'TP53',		1),
	('%s/reference_files/reannoted/Lymphome_B_IAD119887_231_Designed.with_NM.bed' % pipeline_folder	,'Lymphome_B',	1),
	('%s/reference_files/reannoted/Lymphome_T_IAD120574_238_Designed.with_NM.bed' % pipeline_folder	,'Lymphome_T',	1),
	('%s/reference_files/reannoted/IAD83112_241_Designed.with_NM.bed' % pipeline_folder				,'Leuc',		1),
	('%s/reference_files/reannoted/ABL1_NM_005157_Designed.with_NM.bed' % pipeline_folder			,'ABL1',		1),
	('%s/reference_files/reannoted/FLT3_IAD161204_182_Designed.with_NM.bed' % pipeline_folder		,'FLT3',		1),
	('%s/reference_files/reannoted/IAD78219_237_Designed_with_NM.bed' % pipeline_folder				,'LAM',			1),
	('%s/reference_files/reannoted/IAD37093_Designed_with_NM.bed' % pipeline_folder					,'LAM',			2),
	('%s/reference_files/reannoted/IAD62716_182_Designed_with_NM.bed' % pipeline_folder				,'LAM',			3),
	('%s/reference_files/reannoted/IAD165023_231_Designed.with_NM.bed' % pipeline_folder			,'SBT',			9),
	('%s/reference_files/reannoted/IAD154118_231_Designed.with_NM.bed' % pipeline_folder			,'SBT',			8),
	('%s/reference_files/reannoted/IAD119108_231_Designed.with_NM.bed' % pipeline_folder			,'SBT',			7),
	('%s/reference_files/reannoted/IAD108862_231_Designed.with_NM.bed' % pipeline_folder			,'SBT',			5),
	('%s/reference_files/reannoted/IAD94971_233_Designed.with_NM.bed' % pipeline_folder				,'SBT',			4),
	('%s/reference_files/reannoted/IAD72953_231_Designed.with_NM.bed' % pipeline_folder				,'SBT',			3),
	('%s/reference_files/reannoted/SureSelect-HEMATO-v5.sorted.annotated.bed' % pipeline_folder		,'TEST',		1),
	('%s/reference_files/Target-Myeloid_v1-SureSelect.roi.anno.bed' % pipeline_folder				,'LAM-illumina',1),
	('%s/reference_files/Target-Myeloid_capture_v2.roi.anno.padding5.bed' % pipeline_folder			,'LAM-illumina',2),
]

db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

for p in panels:
	panel_path = p[0]
	panel_project = p[1]
	panel_version = p[2]
	panel_name = panel_path.split('/')[-1]#.replace('.annotated.bed','.bed')
	db_cur.execute("SELECT * FROM Panel WHERE panelID='%s'" % panel_name)
	if db_cur.fetchone():
		print "- %s already in database" % panel_name
		continue
	#panel_size = subprocess.check_output("cat %s | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'" % panel_path, shell=True).strip()
	#db_cur.execute("INSERT INTO Panel (panelID, panelProject, panelVersion, panelSize) VALUES ('%s', '%s', %s, %s)" % (panel_name,panel_project,panel_version,panel_size))
	db_cur.execute("INSERT INTO Panel (panelID, panelProject) VALUES ('%s', '%s')" % (panel_name,panel_project))

	panel_reader = open(panel_path,'r')
	chr_count = {'chr1':[],'chr2':[],'chr3':[],'chr4':[],'chr5':[],'chr6':[],'chr7':[],'chr8':[],'chr9':[],'chr10':[],'chr11':[],'chr12':[],'chr13':[],'chr14':[],'chr15':[],'chr16':[],'chr17':[],'chr18':[],'chr19':[],'chr20':[],'chr21':[],'chr22':[],'chrX':[],'chrY':[]}
	base_count = 0
	for line in panel_reader:
		if line.startswith('track'):
			continue
		print line
		line = line.replace('\n','').split('\t')
		chromosome = line[0]
		if 'ABL1' in chromosome:
			chromosome = 'chr9'
		start = int(line[1])
		stop = int(line[2])
		for position in range(start+1,stop+1):
			if position not in chr_count[chromosome]:
				chr_count[chromosome].append(position)
				base_count += 1
		targetedregion_name = line[3]
		random_uuid = uuid.uuid1()
		targetedregion_id = 'T-'+random_uuid.hex[:8]
		gene = line[7].split('GENE=')[-1].split(';')[0]
		details = line[7].split('DETAILS=')[-1].split(';')[0]
		transcript = line[7].split('TRANSCRIPT=')[-1].split(';')[0]
		strand = line[7].split('STRAND=')[-1].split(';')[0]
		nc = NC_data[chromosome]
		db_cur.execute("INSERT INTO TargetedRegion (targetedRegionID, targetedRegionName, panel, chromosome, start, stop, gene, details) VALUES ('%s', '%s', '%s', '%s', %s, %s, '%s', '%s')" % (targetedregion_id,targetedregion_name,panel_name,chromosome,start,stop,gene,details))
		
		# SI GENE ABSENT de la BDD, AJOUTER
		db_cur.execute("SELECT * FROM Gene WHERE geneID='%s'" % gene)
		db_gene = db_cur.fetchone()
		if not db_gene:
			refGene_file.seek(0)
			for rgline in refGene_reader:
				if rgline[1].split('.')[0] == transcript and ('_' not in rgline[2]) :
					print "\t OK"
					transcriptionStart = rgline[4]
					transcriptionStop = rgline[5]	# rgline[6] et rgline[7] sont CodingRegionStart et codingRegionStop
					exons = int(rgline[8])
					exonsStart = rgline[9]
					exonsStop =  rgline[10]
					break
			
			gene_tx = hdp.get_tx_for_gene(gene)
			version = 0
			for item in gene_tx:
				if (item[3].split('.')[0] == transcript) and (item[4] == nc):
					version = max(version,int(item[3].split('.')[-1]))	
			db_cur.execute("INSERT INTO Gene (geneID, chromosome, transcriptionStart, transcriptionStop, strand, exons, exonsStart, exonsStop, transcript, transcriptVersion, NC) VALUES ('%s', '%s', %s, %s, '%s', %s, '%s', '%s', '%s', %s, '%s')" % (gene,chromosome,transcriptionStart,transcriptionStop,strand,exons,exonsStart,exonsStop,transcript,version,nc))
		else:
			# VERIFIER TRANSCRIPT
			if transcript != db_gene['transcript']:
				print "\t-warning : %s already in DB, but transcript different (%s in DB, %s in bed)" % (gene,db_gene['transcript'],transcript)
		
	db_cur.execute("UPDATE Panel SET panelSize=%s WHERE panelID='%s'" % (base_count,panel_name))
	print "- adding : %s, %s, %s b" % (panel_path, panel_project, base_count)
		
db_con.commit()
db_con.close()
