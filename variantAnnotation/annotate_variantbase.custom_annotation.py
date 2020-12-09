import sqlite3
import json
import time
import xlrd
import uuid
import os
import csv
import sys

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

custom_false_positives = global_param['False_positives']
custom_drug_sensitivity = global_param['SBT_sensitivity']
custom_highlight = global_param['LAM_hotspot_highlight']
custom_tp53_umd = global_param['TP53_UMD_variants_EU']
custom_lymphome_controls = global_param['LymphomeT_controls']

## FAIRE OPTIONS? : SI --FULL ON SUPPRIME TOUTES LES ENTREEES POUR LES RECREES ET FAIRE DES "INSERT"
## SI --NEW, FAIRE DES "INSERT OR IGNORE"
print "- Cleaning all VariantAnnotation entries..."
db_cur.execute("DELETE FROM userComment")

#############################
# ADDING CUSTOM ANNOTATIONS #
#############################

print "\n- [%s] adding custom annotations ..." % (time.strftime("%H:%M:%S"))

db_cur.execute("SELECT DISTINCT panelID FROM Panel INNER JOIN Analysis ON Analysis.panel = Panel.panelID INNER JOIN Run ON Run.runID = Analysis.run WHERE platform='ion torrent'")
ion_torrent_panels = ','.join([str(item['panelID']) for item in db_cur.fetchall()])
db_cur.execute("SELECT DISTINCT panelID FROM Panel INNER JOIN Analysis ON Analysis.panel = Panel.panelID INNER JOIN Run ON Run.runID = Analysis.run WHERE platform='illumina'")
illumina_panels = ','.join([str(item['panelID']) for item in db_cur.fetchall()])

# KNOWN FALSE-POSITIVES
print "\t- [%s] adding FalsePositives ..." % (time.strftime("%H:%M:%S"))
with open(custom_false_positives,'r') as fp:
	fp_reader = csv.reader(fp,delimiter='\t')
	for fp_line in fp_reader :
		transcript_without_version = fp_line[0]
		cpos = fp_line[1]
		comment = fp_line[2].replace("'"," ")
		panels = fp_line[3].replace('ion-torrent-panels',ion_torrent_panels).replace('illumina-panels',illumina_panels)
		panels = panels.split(',')

		db_cur.execute("SELECT variantAnnotationID FROM VariantAnnotation WHERE transcript LIKE '%s.%%' AND (transcriptDescription='%s' OR annovarTranscriptDescription='%s')" % (transcript_without_version,cpos,cpos))
		va_ids = db_cur.fetchall()
		for va_id in va_ids:
			va_id = va_id['variantAnnotationID']
			for panel in panels:
				random_uuid = uuid.uuid1()
				usercomment_id = 'C-'+random_uuid.hex[:8]
				db_cur.execute("INSERT INTO UserComment (userCommentID,variantAnnotation,panel,userComment) VALUES ('%s','%s','%s','%s')" % (usercomment_id,va_id,panel,comment))

# DRUG SENSITIVITY
print "\t- [%s] adding Drug Sensitivity ..." % (time.strftime("%H:%M:%S"))
sensitivity_xls = xlrd.open_workbook(custom_drug_sensitivity)
sensitivity_sheet = sensitivity_xls.sheet_by_index(0)
for row in range(sensitivity_sheet.nrows-1):
	transcript_without_version = str(sensitivity_sheet.cell(row+1,4).value)
	cpos = str(sensitivity_sheet.cell(row+1,7).value)
	ppos = str(sensitivity_sheet.cell(row+1,8).value)
	actionability = str(sensitivity_sheet.cell(row+1,9).value)
	if ppos not in ['','p.(=)','p.?','p.=']:
		db_cur.execute("SELECT variantAnnotationID FROM VariantAnnotation WHERE transcript LIKE '%s.%%' AND (proteinDescription='%s' OR annovarProteinDescription='%s')" % (transcript_without_version,ppos,ppos))
	elif cpos != '':
		db_cur.execute("SELECT variantAnnotationID FROM VariantAnnotation WHERE transcript LIKE '%s.%%' AND (transcriptDescription='%s' OR annovarTranscriptDescription='%s')" % (transcript_without_version,cpos,cpos))
	va_ids = db_cur.fetchall()
	for va_id in va_ids:
		variantAnnotationID = va_id['variantAnnotationID']
		db_cur.execute("UPDATE VariantAnnotation SET actionability='%s' WHERE variantAnnotationID='%s'" % (actionability,variantAnnotationID))

## cas particulier : splicing MET
print "\t- [%s] adding Splice MET ..." % (time.strftime("%H:%M:%S"))
db_cur.execute("""SELECT variantAnnotationID,actionability FROM VariantAnnotation 
INNER JOIN Transcript ON Transcript.transcriptID = VariantAnnotation.transcript 
INNER JOIN Variant ON Variant.variantID = VariantAnnotation.variant
WHERE gene='MET' AND region='intronic' AND ((genomicStart BETWEEN 116411873 AND 116411902) OR (genomicStart BETWEEN 116412044 AND 116412087))""")
va_ids = db_cur.fetchall()
for va_id in va_ids:
	variantAnnotationID = va_id['variantAnnotationID']
	db_cur.execute("UPDATE VariantAnnotation SET exon=14 WHERE variantAnnotationID='%s'" % variantAnnotationID)
	if va_id['actionability'] == None :
		db_cur.execute("UPDATE VariantAnnotation SET actionability='MET intron 13/14' WHERE variantAnnotationID='%s'" % variantAnnotationID)


# LAM HOTSPOT HIGHLIGHTS
print "\t- [%s] adding LAM hotspot Highlights ..." % (time.strftime("%H:%M:%S"))
with open(custom_highlight,'r') as hl:
	hl_reader = csv.reader(hl,delimiter='\t')
	for hl_row in hl_reader:
		transcript_without_version = hl_row[0]
		cpos = hl_row[1]

		db_cur.execute("SELECT variantAnnotationID FROM VariantAnnotation WHERE transcript LIKE '%s.%%' AND (transcriptDescription='%s' OR annovarTranscriptDescription='%s')" % (transcript_without_version,cpos,cpos))
		va_ids = db_cur.fetchall()
		for va_id in va_ids:
			va_id = va_id['variantAnnotationID']
			db_cur.execute("UPDATE VariantAnnotation SET highlight=1 WHERE variantAnnotationID='%s'" % va_id)


# TP53 UMD
print "\t- [%s] adding TP53 UMD ..." % (time.strftime("%H:%M:%S"))
with open(custom_tp53_umd,'r') as tp53_umd:
	tp53_umd_reader = csv.reader(tp53_umd,delimiter='\t')
	for tp53_row in tp53_umd_reader:
		cpos = tp53_row[0]
		ppos = tp53_row[1]
		pathoUMD = tp53_row[2]
		commentUMD = tp53_row[3]
		if ppos not in ['','p.(=)','p.?','p.=']:
			db_cur.execute("""SELECT variantAnnotationID FROM VariantAnnotation 
			INNER JOIN Transcript ON Transcript.transcriptID = VariantAnnotation.transcript
			WHERE gene='TP53' AND (proteinDescription='%s' OR annovarProteinDescription='%s')""" % (ppos,ppos))
		elif cpos not in ['','c.(=)','c.?','c.=']:
			db_cur.execute("""SELECT variantAnnotationID FROM VariantAnnotation 
			INNER JOIN Transcript ON Transcript.transcriptID = VariantAnnotation.transcript
			WHERE gene='TP53' AND (transcriptDescription='%s' OR annovarTranscriptDescription='%s')""" % (cpos,cpos))
		va_ids = db_cur.fetchall()
		for va_id in va_ids:
			variantAnnotationID = va_id['variantAnnotationID']
			db_cur.execute("UPDATE VariantAnnotation SET pathoUMD='%s',commentUMD='%s' WHERE variantAnnotationID='%s'" % (pathoUMD,commentUMD,variantAnnotationID))


# LYMPHOME EFS CONTROLS
print "\t- [%s] adding Lymphome EFS Controls ..." % (time.strftime("%H:%M:%S"))
with open(custom_lymphome_controls,'r') as lc:
	lc_reader = csv.reader(lc,delimiter='\t')
	for lc_row in lc_reader :
		transcript_without_version = lc_row[0]
		cpos = lc_row[1]
		nb = len(lc_row[2].split(','))
		comment = "Found in normal dna (x%s)" % nb # "Found in %s controls (%s)" % (nb,lc_row[2])

		db_cur.execute("SELECT variantAnnotationID FROM VariantAnnotation WHERE transcript LIKE '%s.%%' AND (transcriptDescription='%s' OR annovarTranscriptDescription='%s')" % (transcript_without_version,cpos,cpos))
		va_ids = db_cur.fetchall()
		for va_id in va_ids:
			va_id = va_id['variantAnnotationID']
			for panel in ['Lymphome_B','Lymphome_T']:
				random_uuid = uuid.uuid1()
				usercomment_id = 'C-'+random_uuid.hex[:8]
				db_cur.execute("INSERT INTO UserComment (userCommentID,variantAnnotation,panel,userComment) VALUES ('%s','%s','%s','%s')" % (usercomment_id,va_id,panel,comment))

db_con.commit()
db_con.close()
