#!/usr/bin/python
import os
import sys
import csv
import json
import sqlite3
import subprocess

def replace_empty_with_null(to_db):
	mod_to_db = []
	for entry in to_db:
		mod_entry = []
		for column_value in entry:
			if column_value == '':
				mod_entry.append(None)
			else:
				mod_entry.append(column_value)
		mod_to_db.append(tuple(mod_entry))
	return mod_to_db

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

 # __   __   ___      ___  ___          ___          __   __  
# /  ` |__) |__   /\   |  |__     |\ | |__  |  |    |  \ |__) 
# \__, |  \ |___ /~~\  |  |___    | \| |___ |/\|    |__/ |__) 
                                                            
print "- create new db"
new_db_path = global_param['VariantBase'].replace('.db','.new.db')
if os.path.exists(new_db_path):
	subprocess.call(['rm',new_db_path])

db_con = sqlite3.connect(new_db_path)
db_cur = db_con.cursor()

create_run_table = """
CREATE TABLE `Run` (
	`runID`	TEXT NOT NULL,
	`platform`	TEXT NOT NULL,
	`system`	TEXT NOT NULL,
	`runPath`	TEXT,
	`runDate`	TEXT NOT NULL,
	PRIMARY KEY(`runID`)
);"""

create_sample_table = """
CREATE TABLE `Sample` (
	`sampleID`	TEXT NOT NULL,
	`sampleName`	TEXT NOT NULL,
	`gender`	TEXT,
	`pathology`	TEXT,
	`isControl`	INTEGER NOT NULL DEFAULT 0,
	`lims`	TEXT,
	PRIMARY KEY(`sampleID`)
);"""

create_panel_table = """
CREATE TABLE `Panel` (
	`panelID`	TEXT NOT NULL,
	`bedName`	TEXT NOT NULL,
	`panelProject`	TEXT NOT NULL,
	`panelSubProject`	TEXT,
	`panelSize`	INTEGER,
	PRIMARY KEY(`panelID`)
);"""

create_transcript_table = """
CREATE TABLE `Transcript` (
	`transcriptID`	TEXT NOT NULL,
	`gene`	TEXT,
	`chromosome`	TEXT,
	`transcriptionStart`	INTEGER,
	`transcriptionStop`	NUMERIC,
	`exons`	INTEGER,
	`exonsStart`	TEXT,
	`exonsStop`	TEXT,
	`strand`	TEXT,
	`NC`	TEXT,
	PRIMARY KEY(`transcriptID`)
);"""

create_targetedregion_table = """
CREATE TABLE `TargetedRegion` (
	`targetedRegionID`	TEXT NOT NULL,
	`targetedRegionName`	TEXT NOT NULL,
	`chromosome`	TEXT NOT NULL,
	`start`	INTEGER NOT NULL,
	`stop`	INTEGER NOT NULL,
	`details`	TEXT,
	`pool`	INTEGER,
	`panel`	TEXT NOT NULL,
	`transcript`	TEXT NOT NULL,
	PRIMARY KEY(`targetedRegionID`),
	FOREIGN KEY(`panel`) REFERENCES `Panel`(`panelID`),
	FOREIGN KEY(`transcript`) REFERENCES `Transcript`(`transcriptID`)
);"""

create_variant_table = """
CREATE TABLE `Variant` (
	`variantID`	TEXT NOT NULL,
	`chromosome`	TEXT NOT NULL,
	`genomicStart`	INTEGER NOT NULL,
	`genomicStop`	INTEGER NOT NULL,
	`referenceAllele`	TEXT NOT NULL,
	`alternativeAllele`	TEXT NOT NULL,
	`genomicDescription`	TEXT,
	PRIMARY KEY(`variantID`)
);"""

create_variantannotation_table = """
CREATE TABLE `VariantAnnotation` (
	`variantAnnotationID`	TEXT NOT NULL,
	`variant`	TEXT NOT NULL,
	`transcript`	TEXT,
	`commentaire`	TEXT,
	`transcriptDescription`	TEXT,
	`proteinDescription`	TEXT,
	`region`	TEXT,
	`consequence`	TEXT,
	`exon`	INTEGER,
	`intron`	INTEGER,
	`variantType`	TEXT,
	`actionability`	TEXT,
	`intervar`	TEXT,
	`clinvar`	TEXT,
	`cosmic`	TEXT,
	`dbsnp`	TEXT,
	`gnomad`	TEXT,
	`milleGall`	TEXT,
	`milleGeur`	TEXT,
	`nci60`	TEXT,
	`esp`	TEXT,
	`exac`	TEXT,
	`sift`	TEXT,
	`polyphen2`	TEXT,
	`provean`	TEXT,
	`pubmed`	TEXT,
	`vep_consequence`	TEXT,
	`vep_impact`	TEXT,
	`vep_diff`	TEXT,
	`highlight`	INTEGER,
	`pathoUMD`	TEXT,
	`commentUMD`	TEXT,
	`annovarTranscriptDescription`	TEXT,
	`annovarProteinDescription`	TEXT,
	`hgvs`	TEXT,
	`hgvsInfo`	TEXT,
	`annoWarning`	TEXT,
	`lastUpdate`	TEXT,
	PRIMARY KEY(`variantAnnotationID`),
	FOREIGN KEY(`variant`) REFERENCES `Variant`(`variantID`),
	FOREIGN KEY(`transcript`) REFERENCES `Transcript`(`transcriptID`)
);"""

create_analysis_table = """
CREATE TABLE `Analysis` (
	`analysisID`	TEXT NOT NULL,
	`sample`	TEXT NOT NULL,
	`barcode`	TEXT NOT NULL,
	`run`	TEXT NOT NULL,
	`panel`	TEXT NOT NULL,
	`bamPath`	TEXT,
	`analysisDate`	TEXT,
	PRIMARY KEY(`analysisID`),
	FOREIGN KEY(`run`) REFERENCES `Run`(`runID`),
	FOREIGN KEY(`panel`) REFERENCES `Panel`(`panelID`),
	FOREIGN KEY(`sample`) REFERENCES `Sample`(`sampleID`)
);"""

create_variantmetrics_table = """
CREATE TABLE `VariantMetrics` (
	`variantMetricsID`	TEXT NOT NULL,
	`variant`	TEXT NOT NULL,
	`analysis`	TEXT NOT NULL,
	`positionReadDepth`	INTEGER NOT NULL,
	`variantReadDepth`	INTEGER NOT NULL,
	`variantCallingTool`	TEXT NOT NULL,
	`call`	TEXT,
	FOREIGN KEY(`variant`) REFERENCES `Variant`(`variantID`),
	PRIMARY KEY(`variantMetricsID`),
	FOREIGN KEY(`analysis`) REFERENCES `Analysis`(`analysisID`)
);"""

create_usercomment_table = """
CREATE TABLE `UserComment` (
	`userCommentID`	TEXT NOT NULL,
	`variantAnnotation`	TEXT NOT NULL,
	`panel`	TEXT NOT NULL,
	`userComment`	TEXT,
	PRIMARY KEY(`userCommentID`),
	FOREIGN KEY(`variantAnnotation`) REFERENCES `VariantAnnotation`(`variantAnnotationID`),
	FOREIGN KEY(`panel`) REFERENCES `Panel`(`panelID`)
);"""

db_cur.execute(create_run_table)
db_cur.execute(create_sample_table)
db_cur.execute(create_panel_table)
db_cur.execute(create_transcript_table)
db_cur.execute(create_targetedregion_table)
db_cur.execute(create_variant_table)
db_cur.execute(create_variantannotation_table)
db_cur.execute(create_analysis_table)
db_cur.execute(create_variantmetrics_table)
db_cur.execute(create_usercomment_table)

 # ___      __   __   __  ___     __        __      __   __  
# |__  \_/ |__) /  \ |__)  |     /  \ |    |  \    |  \ |__) 
# |___ / \ |    \__/ |  \  |     \__/ |___ |__/    |__/ |__) 
                                                           
print "- export old db"
# first : 'Run','Sample','Panel','Transcript','TargetedRegion','Variant'
# second : 'VariantAnnotation','Analysis','VariantMetrics'
# third : userComment

table_list = ['Run','Sample','Panel','Transcript','TargetedRegion','Variant','VariantAnnotation','Analysis','VariantMetrics','UserComment']
for table in table_list:
	print "\t- export %s" % table
	sqlrequest = "select * from %s;" % table
	subprocess.call(['sqlite3','-header','-csv',global_param['VariantBase'],sqlrequest],stdout=open('/tmp/%s.csv' % table,'w'))

# # DOS2UNIX
# for csvfile in ['Run_old','Run_new','Sample_old','Sample_new','Panel','Transcript','TargetedRegion','Variant_old','Variant_new','VariantAnnotation_old','VariantAnnotation_new','Analysis_old','Analysis_new','VariantMetrics_old','VariantMetrics_new']:
	# subprocess.call(['python','%s/scripts/dos2unix.py' % pipeline_folder,'/tmp/%s.csv' % csvfile,'/tmp/%s.csv' % csvfile])

#          __   __   __  ___     __        __      __   __     ___  __           ___          __   __  
# |  |\/| |__) /  \ |__)  |     /  \ |    |  \    |  \ |__)     |  /  \    |\ | |__  |  |    |  \ |__) 
# |  |  | |    \__/ |  \  |     \__/ |___ |__/    |__/ |__)     |  \__/    | \| |___ |/\|    |__/ |__) 

print "- import old db to new db"
# RUN
with open('/tmp/Run.csv','r') as fin:
	print "\t- import Run"
	dr = csv.DictReader(fin)
	to_db = [(i['runID'], i['platform'], i['system'], i['runPath'], i['runDate']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Run (runID, platform, system, runPath, runDate) VALUES (?, ?, ?, ?, ?);", to_db)
# SAMPLE
with open('/tmp/Sample.csv','r') as fin:
	print "\t- import Sample"
	dr = csv.DictReader(fin)
	to_db = [(i['sampleID'], i['sampleName'], i['gender'], i['pathology'], i['isControl']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Sample (sampleID, sampleName, gender, pathology, isControl) VALUES (?, ?, ?, ?, ?);", to_db)
# PANEL
with open('/tmp/Panel.csv','r') as fin:
	print "\t- import Panel"
	dr = csv.DictReader(fin)
	to_db = [(i['panelID'], i['bedName'], i['panelProject'], i['panelSubProject'], i['panelSize']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Panel (panelID, bedName, panelProject, panelSubProject, panelSize) VALUES (?, ?, ?, ?, ?);", to_db)
# TRANSCRIPT
with open('/tmp/Transcript.csv','r') as fin:
	print "\t- import Transcript"
	dr = csv.DictReader(fin)
	to_db = [(i['transcriptID'], i['gene'], i['chromosome'], i['transcriptionStart'], i['transcriptionStop'], i['strand'], i['exons'], i['exonsStart'], i['exonsStop'], i['NC']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Transcript (transcriptID, gene, chromosome, transcriptionStart, transcriptionStop, strand, exons, exonsStart, exonsStop, NC) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", to_db)
# TARGETED REGION
with open('/tmp/TargetedRegion.csv','r') as fin:
	print "\t- import TargetedRegion"
	dr = csv.DictReader(fin)
	to_db = [(i['targetedRegionID'], i['targetedRegionName'], i['panel'], i['chromosome'], i['start'], i['stop'], i['transcript'], i['details'], i['pool']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO TargetedRegion (targetedRegionID, targetedRegionName, panel, chromosome, start, stop, transcript, details, pool) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);", to_db)
# VARIANT
with open('/tmp/Variant.csv','r') as fin:
	print "\t- import Variant"
	dr = csv.DictReader(fin)
	to_db = [(i['variantID'], i['chromosome'], i['genomicStart'], i['genomicStop'], i['referenceAllele'], i['alternativeAllele']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Variant (variantID,chromosome,genomicStart,genomicStop,referenceAllele,alternativeAllele) VALUES (?, ?, ?, ?, ?, ?);", to_db)
# VARIANT ANNOTATION
with open('/tmp/VariantAnnotation.csv','r') as fin:
	print "\t- import VariantAnnotation"
	dr = csv.DictReader(fin)
	to_db = [(i['variantAnnotationID'], i['variant'], i['transcript'], i['commentaire'], i['transcriptDescription'], i['proteinDescription'], i['region'], i['consequence'], i['exon'], i['intron'], i['actionability'], i['intervar'], i['clinvar'], i['cosmic'], i['dbsnp'], i['gnomad'], i['milleGall'], i['milleGeur'], i['nci60'], i['esp'], i['exac'], i['sift'], i['polyphen2'], i['provean'], i['pubmed'], i['vep_consequence'], i['vep_impact'], i['vep_diff'], i['highlight'], i['pathoUMD'], i['commentUMD'], i['annovarTranscriptDescription'], i['annovarProteinDescription'], i['hgvs'], i['hgvsInfo'], i['annoWarning'], i['lastUpdate']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO VariantAnnotation (variantAnnotationID,variant,transcript,commentaire,transcriptDescription,proteinDescription,region,consequence,exon,intron,actionability,intervar,clinvar,cosmic,dbsnp,gnomad,milleGall,milleGeur,nci60,esp,exac,sift,polyphen2,provean,pubmed,vep_consequence,vep_impact,vep_diff,highlight,pathoUMD,commentUMD,annovarTranscriptDescription,annovarProteinDescription,hgvs,hgvsInfo,annoWarning,lastUpdate) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", to_db)
# ANALYSIS
with open('/tmp/Analysis.csv','r') as fin:
	print "\t- import Analysis"
	dr = csv.DictReader(fin)
	to_db = [(i['analysisID'], i['sample'], i['barcode'], i['run'], i['panel'], i['bamPath'], i['analysisDate']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT INTO Analysis (analysisID, sample, barcode, run, panel, bamPath, analysisDate) VALUES (?, ?, ?, ?, ?, ?, ?);", to_db)
# VARIANT METRICS
with open('/tmp/VariantMetrics.csv','r') as fin:
	print "\t- import VariantMetrics"
	dr = csv.DictReader(fin)
	to_db = [(i['variantMetricsID'], i['variant'], i['analysis'], i['positionReadDepth'], i['variantReadDepth'], i['variantCallingTool'], i['call']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO VariantMetrics (variantMetricsID, variant, analysis, positionReadDepth, variantReadDepth, variantCallingTool, call) VALUES (?, ?, ?, ?, ?, ?, ?);", to_db)
# USERCOMMENT
with open('/tmp/UserComment.csv','r') as fin:
	print "\t- import UserComment"
	dr = csv.DictReader(fin)
	to_db = [(i['userCommentID'], i['panel'], i['variantAnnotation'], i['userComment']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO UserComment (userCommentID, panel, variantAnnotation, userComment) VALUES (?, ?, ?, ?);", to_db)

db_con.commit()
db_con.close()

 # __   ___  __             __   ___     __        __      __   __      __               ___      
# |__) |__  |__) |     /\  /  ` |__     /  \ |    |  \    |  \ |__)    |__) \ /    |\ | |__  |  | 
# |  \ |___ |    |___ /~~\ \__, |___    \__/ |___ |__/    |__/ |__)    |__)  |     | \| |___ |/\| 
                                                                                                
print "- replace old db by new"
subprocess.call(['mv',global_param['VariantBase'],global_param['VariantBase'].replace('.db','.old.db')])
subprocess.call(['mv',new_db_path,global_param['VariantBase']])

print "- done."
