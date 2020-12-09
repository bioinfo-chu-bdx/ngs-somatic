#!/usr/bin/python
import os
import sys
import csv
import json
import sqlite3

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

new_db_path = global_param['VariantBase'].replace('.db','.new.db')

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

db_con.commit()
db_con.close()
