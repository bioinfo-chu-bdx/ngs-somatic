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

new_db_path = global_param['VariantBase'].replace('.db','.new.db')

db_con = sqlite3.connect(new_db_path)
db_cur = db_con.cursor()

# DOS2UNIX
for csvfile in ['Run_old','Run_new','Sample_old','Sample_new','Panel','Transcript','TargetedRegion','Variant_old','Variant_new','VariantAnnotation_old','VariantAnnotation_new','Analysis_old','Analysis_new','VariantMetrics_old','VariantMetrics_new']:
	subprocess.call(['python','%s/scripts/dos2unix.py' % pipeline_folder,'/media/stuff/variantBase_export/%s.csv' % csvfile,'/media/stuff/variantBase_export/%s.csv' % csvfile])

# RUN
with open('/media/stuff/variantBase_export/Run_old.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['runID'], i['platform'], i['system'], i['runPath'], i['runDate']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Run (runID, platform, system, runPath, runDate) VALUES (?, ?, ?, ?, ?);", to_db)
with open('/media/stuff/variantBase_export/Run_new.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['runID'], i['platform'], i['system'], i['runPath'], i['runDate']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Run (runID, platform, system, runPath, runDate) VALUES (?, ?, ?, ?, ?);", to_db)

# SAMPLE
with open('/media/stuff/variantBase_export/Sample_old.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['sampleID'], i['sampleName'], i['gender'], i['pathology'], i['isControl']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Sample (sampleID, sampleName, gender, pathology, isControl) VALUES (?, ?, ?, ?, ?);", to_db)
with open('/media/stuff/variantBase_export/Sample_new.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['sampleID'], i['sampleName'], i['gender'], i['pathology'], i['isControl']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Sample (sampleID, sampleName, gender, pathology, isControl) VALUES (?, ?, ?, ?, ?);", to_db)

# PANEL
with open('/media/stuff/variantBase_export/Panel.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['panelID'], i['bedName'], i['panelProject'], i['panelSubProject'], i['panelSize']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Panel (panelID, bedName, panelProject, panelSubProject, panelSize) VALUES (?, ?, ?, ?, ?);", to_db)

# # GENE
# with open('/media/stuff/variantBase_export/Gene.csv','r') as fin:
	# dr = csv.DictReader(fin)
	# to_db = [(i['geneID'], i['chromosome'], i['transcriptionStart'], i['transcriptionStop'], i['strand'], i['exons'], i['exonsStart'], i['exonsStop'], i['transcript'], i['transcriptVersion'], i['NC']) for i in dr]
	# to_db = replace_empty_with_null(to_db)
# db_cur.executemany("INSERT OR IGNORE INTO Gene (geneID, chromosome, transcriptionStart, transcriptionStop, strand, exons, exonsStart, exonsStop, transcript, transcriptVersion, NC) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", to_db)

# TRANSCRIPT
with open('/media/stuff/variantBase_export/Transcript.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['transcriptID'], i['gene'], i['chromosome'], i['transcriptionStart'], i['transcriptionStop'], i['strand'], i['exons'], i['exonsStart'], i['exonsStop'], i['NC']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Transcript (transcriptID, gene, chromosome, transcriptionStart, transcriptionStop, strand, exons, exonsStart, exonsStop, NC) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", to_db)

# TARGETED REGION
with open('/media/stuff/variantBase_export/TargetedRegion.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['targetedRegionID'], i['targetedRegionName'], i['panel'], i['chromosome'], i['start'], i['stop'], i['transcript'], i['details'], i['pool']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO TargetedRegion (targetedRegionID, targetedRegionName, panel, chromosome, start, stop, transcript, details, pool) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);", to_db)

# VARIANT
with open('/media/stuff/variantBase_export/Variant_old.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['variantID'], i['chromosome'], i['genomicStart'], i['genomicStop'], i['referenceAllele'], i['alternativeAllele']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Variant (variantID,chromosome,genomicStart,genomicStop,referenceAllele,alternativeAllele) VALUES (?, ?, ?, ?, ?, ?);", to_db)
with open('/media/stuff/variantBase_export/Variant_new.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['variantID'], i['chromosome'], i['genomicStart'], i['genomicStop'], i['referenceAllele'], i['alternativeAllele']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO Variant (variantID,chromosome,genomicStart,genomicStop,referenceAllele,alternativeAllele) VALUES (?, ?, ?, ?, ?, ?);", to_db)

# VARIANT ANNOTATION
with open('/media/stuff/variantBase_export/VariantAnnotation_old.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['variantAnnotationID'], i['variant'], i['transcript'], i['commentaire'], i['transcriptDescription'], i['proteinDescription'], i['region'], i['consequence'], i['exon'], i['intron'], i['actionability'], i['intervar'], i['clinvar'], i['cosmic'], i['dbsnp'], i['gnomad'], i['milleGall'], i['milleGeur'], i['nci60'], i['esp'], i['exac'], i['sift'], i['polyphen2'], i['provean'], i['pubmed'], i['vep_consequence'], i['vep_impact'], i['vep_diff'], i['highlight'], i['pathoUMD'], i['commentUMD'], i['annovarTranscriptDescription'], i['annovarProteinDescription'], i['hgvs'], i['hgvsInfo'], i['annoWarning'], i['lastUpdate']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO VariantAnnotation (variantAnnotationID,variant,transcript,commentaire,transcriptDescription,proteinDescription,region,consequence,exon,intron,actionability,intervar,clinvar,cosmic,dbsnp,gnomad,milleGall,milleGeur,nci60,esp,exac,sift,polyphen2,provean,pubmed,vep_consequence,vep_impact,vep_diff,highlight,pathoUMD,commentUMD,annovarTranscriptDescription,annovarProteinDescription,hgvs,hgvsInfo,annoWarning,lastUpdate) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", to_db)
with open('/media/stuff/variantBase_export/VariantAnnotation_new.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['variantAnnotationID'], i['variant'], i['transcript'], i['commentaire'], i['transcriptDescription'], i['proteinDescription'], i['region'], i['consequence'], i['exon'], i['intron'], i['actionability'], i['intervar'], i['clinvar'], i['cosmic'], i['dbsnp'], i['gnomad'], i['milleGall'], i['milleGeur'], i['nci60'], i['esp'], i['exac'], i['sift'], i['polyphen2'], i['provean'], i['pubmed'], i['vep_consequence'], i['vep_impact'], i['vep_diff'], i['highlight'], i['pathoUMD'], i['commentUMD'], i['annovarTranscriptDescription'], i['annovarProteinDescription'], i['hgvs'], i['hgvsInfo'], i['annoWarning'], i['lastUpdate']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO VariantAnnotation (variantAnnotationID,variant,transcript,commentaire,transcriptDescription,proteinDescription,region,consequence,exon,intron,actionability,intervar,clinvar,cosmic,dbsnp,gnomad,milleGall,milleGeur,nci60,esp,exac,sift,polyphen2,provean,pubmed,vep_consequence,vep_impact,vep_diff,highlight,pathoUMD,commentUMD,annovarTranscriptDescription,annovarProteinDescription,hgvs,hgvsInfo,annoWarning,lastUpdate) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", to_db)

# ANALYSIS
with open('/media/stuff/variantBase_export/Analysis_old.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['analysisID'], i['sample'], i['barcode'], i['run'], i['panel'], i['bamPath'], i['analysisDate']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT INTO Analysis (analysisID, sample, barcode, run, panel, bamPath, analysisDate) VALUES (?, ?, ?, ?, ?, ?, ?);", to_db)
with open('/media/stuff/variantBase_export/Analysis_new.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['analysisID'], i['sample'], i['barcode'], i['run'], i['panel'], i['bamPath'], i['analysisDate']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT INTO Analysis (analysisID, sample, barcode, run, panel, bamPath, analysisDate) VALUES (?, ?, ?, ?, ?, ?, ?);", to_db)

# VARIANT METRICS
with open('/media/stuff/variantBase_export/VariantMetrics_old.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['variantMetricsID'], i['variant'], i['analysis'], i['positionReadDepth'], i['variantReadDepth'], i['variantCallingTool'], i['call']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO VariantMetrics (variantMetricsID, variant, analysis, positionReadDepth, variantReadDepth, variantCallingTool, call) VALUES (?, ?, ?, ?, ?, ?, ?);", to_db)
with open('/media/stuff/variantBase_export/VariantMetrics_new.csv','r') as fin:
	dr = csv.DictReader(fin)
	to_db = [(i['variantMetricsID'], i['variant'], i['analysis'], i['positionReadDepth'], i['variantReadDepth'], i['variantCallingTool'], i['call']) for i in dr]
	to_db = replace_empty_with_null(to_db)
db_cur.executemany("INSERT OR IGNORE INTO VariantMetrics (variantMetricsID, variant, analysis, positionReadDepth, variantReadDepth, variantCallingTool, call) VALUES (?, ?, ?, ?, ?, ?, ?);", to_db)

db_con.commit()
db_con.close()
