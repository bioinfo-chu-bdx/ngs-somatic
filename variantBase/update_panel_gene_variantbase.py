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
panels = {
	# ('%s/reference_files/reannoted/Target_ColonLung_v10_IAD172906_231.bed' % pipeline_folder		,'SBT',			10),
	# ('%s/reference_files/reannoted/LAM2018_IAD143291_236_Designed_with_NM.bed' % pipeline_folder	,'LAM',			4),
	# ('%s/reference_files/reannoted/TP53.20140108.designed_with_NM.bed' % pipeline_folder			,'TP53',		1),
	# ('%s/reference_files/reannoted/Lymphome_B_IAD119887_231_Designed.with_NM.bed' % pipeline_folder	,'Lymphome_B',	1),
	# ('%s/reference_files/reannoted/Lymphome_T_IAD120574_238_Designed.with_NM.bed' % pipeline_folder	,'Lymphome_T',	1),
	# ('%s/reference_files/reannoted/IAD83112_241_Designed.with_NM.bed' % pipeline_folder				,'Leuc',		1),
	# ('%s/reference_files/reannoted/ABL1_NM_005157_Designed.with_NM.bed' % pipeline_folder			,'ABL1',		1),
	# ('%s/reference_files/reannoted/FLT3_IAD161204_182_Designed.with_NM.bed' % pipeline_folder		,'FLT3',		1),
	# ('%s/reference_files/reannoted/IAD78219_237_Designed_with_NM.bed' % pipeline_folder				,'LAM',			1),
	# ('%s/reference_files/reannoted/IAD37093_Designed_with_NM.bed' % pipeline_folder					,'LAM',			2),
	# ('%s/reference_files/reannoted/IAD62716_182_Designed_with_NM.bed' % pipeline_folder				,'LAM',			3),
	# ('%s/reference_files/reannoted/IAD165023_231_Designed.with_NM.bed' % pipeline_folder			,'SBT',			9),
	# ('%s/reference_files/reannoted/IAD154118_231_Designed.with_NM.bed' % pipeline_folder			,'SBT',			8),
	# ('%s/reference_files/reannoted/IAD119108_231_Designed.with_NM.bed' % pipeline_folder			,'SBT',			7),
	# ('%s/reference_files/reannoted/IAD108862_231_Designed.with_NM.bed' % pipeline_folder			,'SBT',			5),
	# ('%s/reference_files/reannoted/IAD94971_233_Designed.with_NM.bed' % pipeline_folder				,'SBT',			4),
	# ('%s/reference_files/reannoted/IAD72953_231_Designed.with_NM.bed' % pipeline_folder				,'SBT',			3),
	# ('%s/reference_files/reannoted/SureSelect-HEMATO-v5.sorted.annotated.bed' % pipeline_folder		,'TEST',		1),
	# ('%s/reference_files/Target-Myeloid_v1-SureSelect.roi.anno.bed' % pipeline_folder				,'LAM-illumina',1),
	'LAM-illumina-v2':{'path':'%s/reference_files/Target-Myeloid_capture_v2.roi.anno.padding5.bed' % pipeline_folder,'project':'LAM','subproject':'panel-capture'},
	'LAM-illumina-v1':{'path':'%s/reference_files/Target-Myeloid_capture_v1.roi.anno.bed' % pipeline_folder,'project':'LAM','subproject':'panel-capture'}
}

db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

db_cur.execute("SELECT * FROM Transcript")
db_transcripts = db_cur.fetchall()
transcript2version_db = {}
for db_transcript in db_transcripts:
	# TODO : il peut y avoir plusieurs transcripts
	transcript_without_version = db_transcript['transcriptID'].split('.')[0]
	version = int(db_transcript['transcriptID'].split('.')[-1])
	if not transcript_without_version in transcript2version_db:
		transcript2version_db[transcript_without_version] = [version]
	else:
		transcript2version_db[transcript_without_version].append(version)

for panel in panels:
	panelID = panel
	print "- %s" % panelID
	bedName = panels[panelID]['path'].split('/')[-1]
	panelProject = panels[panelID]['project']
	panelSubProject = panels[panelID]['subproject']
	# panelSize = subprocess.check_output("cat %s | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'" % panels[panelID]['path'], shell=True).strip()

	# 1 : AJOUTER / UPDATE PANEL ENTRY IN DATABASE
	db_cur.execute("SELECT * FROM Panel WHERE panelID='%s'" % panelID)
	db_panel = db_cur.fetchone()
	if db_panel is not None:
		print "\t - %s already in database" % panelID
		if db_panel['bedName'] != bedName:
			print "\t - updating bedName : %s -> %s" % (db_panel['bedName'],bedName)
			db_cur.execute("UPDATE Panel SET bedName='%s' WHERE panelID='%s'" % (bedName, panelID))
		if db_panel['panelProject'] != panelProject:
			print "\t - updating panelProject : %s -> %s" % (db_panel['panelProject'],panelProject)
			db_cur.execute("UPDATE Panel SET panelProject='%s' WHERE panelID='%s'" % (panelProject, panelID))
		if db_panel['panelSubProject'] != panelSubProject:
			print "\t - updating panelSubProject : %s -> %s" % (db_panel['panelSubProject'],panelSubProject)
			db_cur.execute("UPDATE Panel SET panelSubProject='%s' WHERE panelID='%s'" % (panelSubProject, panelID))
	else:
		db_cur.execute("INSERT INTO Panel (panelID, bedName, panelProject, panelSubProject) VALUES ('%s', '%s', '%s', '%s')" % (panelID, bedName, panelProject, panelSubProject))

	# PARSE BED
	panel_reader = open(panels[panelID]['path'],'r')
	chr_count = {'chr1':[],'chr2':[],'chr3':[],'chr4':[],'chr5':[],'chr6':[],'chr7':[],'chr8':[],'chr9':[],'chr10':[],'chr11':[],'chr12':[],'chr13':[],'chr14':[],'chr15':[],'chr16':[],'chr17':[],'chr18':[],'chr19':[],'chr20':[],'chr21':[],'chr22':[],'chrX':[],'chrY':[]}
	base_count = 0
	z = 0
	for line in panel_reader:
		if line.startswith('track'):
			continue
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
		gene = line[7].split('GENE=')[-1].split(';')[0]
		details = line[7].split('DETAILS=')[-1].split(';')[0]
		transcript_without_version = line[7].split('TRANSCRIPT=')[-1].split(';')[0]
		strand = line[7].split('STRAND=')[-1].split(';')[0]
		nc = NC_data[chromosome]
		print "\t - %s_%s\t%s\t(%s:%s-%s)" % (gene,details,transcript_without_version,chromosome,start,stop)

		# 2 : CHECK TRANSCRIPT PRESENCE IN BDD, AND VERSION NUMBER

		# TRANSCRIPT VERSION
		gene_tx = hdp.get_tx_for_gene(gene)
		version_in_uta = 1
		for item in gene_tx:
			if (item[3].split('.')[0] == transcript_without_version) and (item[4] == nc):
				version_in_uta = max(version_in_uta,int(item[3].split('.')[-1]))

		# TRANSCRIPT INFO
		transcriptionStart = None
		transcriptionStop = None
		exons = None
		exonsStart = None
		exonsStop =  None
		trfound = False
		refGene_file.seek(0)
		for rgline in refGene_reader:
			if rgline[1].split('.')[0] == transcript_without_version and ('_' not in rgline[2]) :
				transcriptionStart = rgline[4]
				transcriptionStop = rgline[5]	# rgline[6] et rgline[7] sont CodingRegionStart et codingRegionStop
				exons = int(rgline[8])
				exonsStart = rgline[9]
				exonsStop =  rgline[10]
				trfound = True
				break

		if not trfound:
			print "ERROR : transcript not found in RefSeq, update manually transcriptionStart, transcriptionStop, exons, exonsStart, exonsStop..."

		# TRANSCRIPT NOT IN DB
		if not (transcript_without_version in transcript2version_db.keys()):
			print "\t\t - transcript %s not found in DB" % transcript_without_version
			transcript = '%s.%s' % (transcript_without_version,version_in_uta)
			print "\t\t\t - adding %s" % transcript
			db_cur.execute("INSERT INTO Transcript (transcriptID, gene, chromosome, transcriptionStart, transcriptionStop, exons, exonsStart, exonsStop, strand, NC) VALUES ('%s', '%s', '%s', %s, %s, %s, '%s', '%s', '%s', '%s')" % (transcript,gene,chromosome,transcriptionStart,transcriptionStop,exons,exonsStart,exonsStop,strand,nc))
			transcript2version_db[transcript_without_version] = [version_in_uta]

		# TRANSCRIPT VERSION NEED AN UPDATE
		elif version_in_uta not in transcript2version_db[transcript_without_version]:
			print "- NEW version for %s : %s->%s, updating" % (transcript_without_version,transcript2version_db[transcript_without_version],version_in_uta)
			transcript = '%s.%s' % (transcript_without_version,version_in_uta)
			db_cur.execute("INSERT INTO Transcript (transcriptID, gene, chromosome, transcriptionStart, transcriptionStop, exons, exonsStart, exonsStop, strand, NC) VALUES ('%s', '%s', '%s', %s, %s, %s, '%s', '%s', '%s', '%s')" % (transcript,gene,chromosome,transcriptionStart,transcriptionStop,exons,exonsStart,exonsStop,strand,nc))
			transcript2version_db[transcript_without_version].append(version_in_uta)
		else:
			transcript = '%s.%s' % (transcript_without_version,version_in_uta)

		# 3 : AJOUTER / UPDATE  TargetedRegion
		db_cur.execute("SELECT * FROM TargetedRegion WHERE panel='%s' AND chromosome='%s' AND start=%s AND stop=%s" % (panelID,chromosome,start,stop))
		db_targetedregion = db_cur.fetchone()
		if db_targetedregion is not None:
			if db_targetedregion['transcript'] != transcript :
				print "\t\t - TargetedRegion %s (%s) found in DB but different, updating transcript %s -> %s" % (db_targetedregion['targetedRegionID'],db_targetedregion['targetedRegionName'],db_targetedregion['transcript'],transcript)
				db_cur.execute("UPDATE TargetedRegion SET transcript='%s', details='%s' WHERE TargetedRegionID='%s'" % (transcript,details,db_targetedregion['targetedRegionID']))
			if db_targetedregion['targetedRegionName'] != targetedregion_name:
				print "\t\t - TargetedRegion %s (%s) found in DB but different, updating targetedRegionName %s -> %s" % (db_targetedregion['targetedRegionID'],db_targetedregion['targetedRegionName'],db_targetedregion['targetedRegionName'],targetedregion_name)
				db_cur.execute("UPDATE TargetedRegion SET targetedRegionName='%s' WHERE TargetedRegionID='%s'" % (targetedregion_name,db_targetedregion['targetedRegionID']))
		else:
			random_uuid = uuid.uuid1()
			targetedregion_id = 'T-'+random_uuid.hex[:8]
			db_cur.execute("INSERT INTO TargetedRegion (targetedRegionID, panel, transcript, targetedRegionName, chromosome, start, stop, details) VALUES ('%s', '%s', '%s', '%s', '%s', %s, %s, '%s')" % (targetedregion_id,panelID,transcript,targetedregion_name,chromosome,start,stop,details))
			print "\r\t - %s new TargetedRegion added" % z,
			z+=1

	if db_panel is not None:
		if db_panel['panelSize'] != base_count:
			print "\t - updating panelSize : %s -> %s" % (db_panel['panelSize'],base_count)
			db_cur.execute("UPDATE Panel SET panelSize=%s WHERE panelID='%s'" % (base_count,panelID))

db_con.commit()
db_con.close()
