#!/usr/bin/python
from optparse import OptionParser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.variantmapper
import hgvs.exceptions
import hgvs.normalizer
import hgvs.exceptions
import hgvs.validator
import hgvs.parser
import sqlite3
import numpy
import json
import time
import os
import sys

# THIS SCRIPT FILL HGVS RESULTS IN VariantBase.db
# USAGE : python annotate_variantbase.step0.py --option

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

def reverse(seq):
	seq = seq.replace('A','W').replace('T','X').replace('G','Y').replace('C','Z')
	seq = seq.replace('W','T').replace('X','A').replace('Y','C').replace('Z','G')
	seq = seq[::-1]
	return seq

### GATHERING PARAMETERS ############################################################

parser = OptionParser()
parser.add_option('-n', '--new',	help="annotate new variants only", 	dest='new', 	default=False, action='store_true')
parser.add_option('-f', '--full',	help="annotate full database", 		dest='full', 	default=False, action='store_true')
parser.add_option('-c', '--clean',	help="delete all variant annotation entries", 		dest='clean', 	default=False, action='store_true')
# parser.add_option('-u', '--update',	help="update old variants", 		dest='update', 	default=False, action='store_true')
(options, args) = parser.parse_args()

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

with open(global_param['NC_data'], 'r') as ncdata:
	NC_data = json.load(ncdata)
db_path = global_param['VariantBase']

db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
hn = hgvs.normalizer.Normalizer(hdp)
hv = hgvs.validator.Validator(hdp)
am37 = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37')
am37no = hgvs.assemblymapper.AssemblyMapper(hdp, normalize=False, assembly_name='GRCh37') # for c to g conversion when c already hgvs

###################################################################################################################################

if options.clean:
	print "- Cleaning all VariantAnnotation entries..."
	db_cur.execute("DELETE FROM VariantAnnotation")

# CORRECT HGVS FOR NEW VARIANTMETRICS OCCURENCE ADDED
print "- checking HGVS in variantMetrics..."
db_cur.execute("""SELECT variantMetricsID,variantID,hgvsInfo FROM VariantMetrics
INNER JOIN Variant ON Variant.variantID = VariantMetrics.variant
INNER JOIN VariantAnnotation ON VariantAnnotation.variant = Variant.variantID
WHERE hgvs = 'no'""")
db_non_hgvs_vms = db_cur.fetchall()
for db_non_hgvs_vm in db_non_hgvs_vms:
	db_cur.execute("SELECT variant FROM VariantAnnotation WHERE variantAnnotationID=?" , (db_non_hgvs_vm['hgvsInfo'],))
	db_hgvs = db_cur.fetchone()
	print "- [HGVS] : Changing variantMetrics %s variant : %s -> %s" % (db_non_hgvs_vm['variantMetricsID'],db_non_hgvs_vm['variantID'],db_hgvs['variant'])
	db_cur.execute("UPDATE VariantMetrics SET variant = ? WHERE VariantMetricsID = ?" , (db_hgvs['variant'],db_non_hgvs_vm['variantMetricsID']))

# FOR ALL VARIANTS AND FOR ALL CORRESPONDING TRANSCRIPTS, CREATE VARIANT ANNOTATION ENTRIES
print "- adding new variantAnnotation entries if needed..."
db_cur.execute("SELECT * FROM Variant")
db_variants = db_cur.fetchall()
db_cur.execute("SELECT * FROM Transcript")
db_transcripts = db_cur.fetchall()
for db_variant in db_variants:
	for db_transcript in db_transcripts:
		if db_variant['chromosome'] == db_transcript['chromosome'] and ((db_transcript['transcriptionStart']-5000)<db_variant['genomicStart']<(db_transcript['transcriptionStop']+5000)):
			variantAnnotationID = '%s:%s' % (db_transcript['transcriptID'],db_variant['variantID'])
			# IF NOT EXISTANT, INSERT NEW VARIANTANNOTATION
			db_cur.execute("SELECT variantAnnotationID FROM VariantAnnotation WHERE variantAnnotationID=?" , (variantAnnotationID,))
			if not db_cur.fetchone():
				print "- new VariantAnnotation entry %s" % variantAnnotationID
				db_cur.execute("INSERT INTO VariantAnnotation (variantAnnotationID, variant, transcript) VALUES (?,?,?)" , (variantAnnotationID, db_variant['variantID'], db_transcript['transcriptID']))

if options.new:
	db_cur.execute("SELECT * FROM VariantAnnotation WHERE hgvs is NULL")#(hgvs is NULL) OR (hgvs='yes' AND region is NULL)")
elif options.full:
	db_cur.execute("SELECT * FROM VariantAnnotation")
else:
	print "error : select either --new or --full"
	exit()

db_variant_annotations = db_cur.fetchall()
print "- %s new variantAnnotation entries to process" % len(db_variant_annotations)
lastUpdate = time.strftime("%Y%m%d")#'%s' % (time.strftime('%d/%m/%Y'))
z=0
for db_variant_annotation in db_variant_annotations:
	z+=1
	variantAnnotationID = db_variant_annotation['variantAnnotationID']
	variantID = db_variant_annotation['variant']
	db_cur.execute("SELECT * FROM Variant WHERE variantID=?" , (variantID,))
	db_variant = db_cur.fetchone()
	chromosome = db_variant['chromosome']
	start = db_variant['genomicStart']
	stop = db_variant['genomicStop']
	ref = db_variant['referenceAllele']
	alt = db_variant['alternativeAllele']
	transcript = db_variant_annotation['transcript']
	db_cur.execute("SELECT * FROM Transcript WHERE transcriptID=?" , (transcript,))
	db_transcript = db_cur.fetchone()
	gene = db_transcript['gene']
	nc = db_transcript['NC']
	strand = db_transcript['strand']

	# GENOMIC DESCRIPTION
	if ref == '-':
		variant_type = 'INS'
		genomicDescription = 'g.%s_%sins%s' % (start,stop,alt)
	elif alt == '-': 
		variant_type = 'DEL'
		if len(ref) > 1:
			genomicDescription = 'g.%s_%sdel%s' % (start,stop,ref)
		else:
			genomicDescription = 'g.%sdel%s' % (start,ref)
	elif len(ref) > 1 or len(alt) > 1:
		variant_type = 'DELINS'
		if len(ref) > 1:
			genomicDescription = 'g.%s_%sdelins%s' % (start,stop,alt)
		else:
			genomicDescription = 'g.%sdelins%s' % (start,alt)
	else: # SNP
		variant_type = 'SNP'
		genomicDescription = 'g.%s%s>%s' % (start,ref,alt)
	completeGenomicDescription = '%s:%s' % (chromosome,genomicDescription)
	db_cur.execute("UPDATE Variant SET genomicDescription=? WHERE variantID=?" , (completeGenomicDescription,variantID))
	db_cur.execute("UPDATE VariantAnnotation SET variantType=? WHERE variantAnnotationID=?" , (variant_type,variantAnnotationID))

	# STEP 1 : PARSE VARIANT (g)
	g0 = hp.parse_hgvs_variant('%s:%s' % (nc,genomicDescription))
	# print "- %s: %s (%s:%s-%s:%s>%s)" % (z,g0,chromosome,start,stop,ref,alt)
	print "- %s: %s" % (z,completeGenomicDescription)

	# STEP 2 : VALIDATE (g)
	try:
		hv.validate(g0)
	except Exception as e:
	# except hgvs.exceptions.HGVSError as e:
		error_message = 'error : %s' % str(e)
		print "\t- %s" % error_message
		db_cur.execute("UPDATE VariantAnnotation SET hgvs='error', hgvsInfo=? WHERE variantAnnotationID=?" , (error_message,variantAnnotationID))
		continue

	# STEP 3 : g to c (with normalization)
	try:
		c = am37.g_to_c(g0,transcript)
	except Exception as e:
	# except hgvs.exceptions.HGVSError as e:
		error_message = 'error : %s' % str(e)
		print "\t- %s" % error_message
		db_cur.execute("UPDATE VariantAnnotation SET hgvs='error', hgvsInfo=? WHERE variantAnnotationID=?" , (error_message,variantAnnotationID))
		continue

	g = am37no.c_to_g(c) # Get coherent g (if not already) from normalyzed c (3' rule)
	# (why dont I normalize g in first place? -> because normalize(g) is useless, it just right-shift in the genomic sens (no 3' rule)

	# STEP 4 : DETERMINATE HGVS. IF HGVS COORDINATES ARE DIFFERENT, CREATE NEW ENTRIES FOR VARIANT AND VARIANTANNOTATION
	if g != g0:
		start = g.posedit.pos.start.base
		end = g.posedit.pos.end.base
		if 'dup' in str(g):
			if strand == 'forward':
				start = g.posedit.pos.end.base
				end = g.posedit.pos.end.base+1
			else:
				start = g.posedit.pos.start.base-1
				end = g.posedit.pos.start.base
			ref = '-'
			alt = str(g.posedit.edit.ref)
			variantType = 'DUP'
		elif 'inv' in str(g):
			ref = str(g.posedit.edit.ref)
			alt = reverse(str(g.posedit.edit.ref))
			variantType = 'INV'
		else:
			ref = str(g.posedit.edit.ref).replace('None','-')
			alt = str(g.posedit.edit.alt).replace('None','-')
			if ref == '-':
				variantType = 'INS'
			elif alt == '-':
				variantType = 'DEL'
			elif len(ref) > 1 or len(alt) > 1: 
				variantType = 'DELINS'
			else:
				variantType = 'SNV'
		print "\t- HGVS should be : %s (%s:%s-%s:%s>%s)" % (g,chromosome,start,end,ref,alt)
		if (start == g0.posedit.pos.start.base) and (end == g0.posedit.pos.end.base):
			genomicDescription = '%s:%s' % (chromosome,str(g).split(':')[-1])
			db_cur.execute("UPDATE Variant SET genomicDescription=? WHERE variantID=? " , (genomicDescription,variantID))
			db_cur.execute("UPDATE VariantAnnotation SET hgvs='yes', variantType=?, lastUpdate=? WHERE variantAnnotationID=?" , (variantType,lastUpdate,variantAnnotationID))

		else: # DIFFERENT COORDINATES FOR HGVS : CREATION OF NEW ENTRIES
			bad_variantID = variantID
			bad_variantAnnotationID = variantAnnotationID
			variantID = '%s:%s-%s:%s>%s' % (chromosome,start,end,ref,alt)
			variantAnnotationID = '%s:%s' % (transcript,variantID)
			genomicDescription = '%s:%s' % (chromosome,str(g).split(':')[-1])

			# NEW VARIANT ENTRY
			db_cur.execute("SELECT * FROM Variant WHERE variantID=?" , (variantID,))
			if db_cur.fetchone() is None:
				print "\t\t- creating new Variant entry in db : %s ..." % variantID
				db_cur.execute("INSERT INTO Variant (variantID, chromosome, genomicStart, genomicStop, referenceAllele, alternativeAllele, genomicDescription) VALUES (?,?,?,?,?,?,?)" , (variantID,chromosome,start,end,ref,alt,genomicDescription))

			# NEW VARIANTANNOTATION ENTRY
			db_cur.execute("SELECT * FROM VariantAnnotation WHERE variantAnnotationID=?" , (variantAnnotationID,))
			if db_cur.fetchone() is None:
				print "\t\t- creating new VariantAnnotation entry in db : %s ..." % variantAnnotationID
				db_cur.execute("INSERT INTO VariantAnnotation (variantAnnotationID, variant, transcript, variantType) VALUES (?,?,?,?)" , (variantAnnotationID,variantID,transcript,variantType))
			db_cur.execute("UPDATE VariantAnnotation SET hgvs='no', hgvsInfo=?, lastUpdate=? WHERE variantAnnotationID = ?" , (variantAnnotationID,lastUpdate,bad_variantAnnotationID))

			# MODIFY VARIANTMETRICS
			db_cur.execute("SELECT * FROM VariantMetrics WHERE variant = ?" , (bad_variantID,))
			varmetrics = db_cur.fetchall()
			for varmetric in varmetrics:
				VariantMetricsID = varmetric['variantMetricsID']
				print "\t\t- updating variantMetrics %s : %s -> %s ..." % (VariantMetricsID,bad_variantID,variantID)
				db_cur.execute("UPDATE VariantMetrics SET variant = ? WHERE VariantMetricsID = ?" , (variantID,VariantMetricsID))

	# STEP 5 : VALIDATE(c)
	c_pos = str(c).split(':')[-1]
	print "\t- transcript description : %s" % str(c)
	db_cur.execute("UPDATE VariantAnnotation SET transcriptDescription=?, hgvs='yes' WHERE variantAnnotationID=?" , (c_pos,variantAnnotationID))
	try:
		hv.validate(c)
	except Exception as e:
	# except hgvs.exceptions.HGVSError as e:
		warning_message = 'warning : %s' % str(e)
		print "\t- %s" % warning_message
		db_cur.execute("UPDATE VariantAnnotation SET hgvs='warning', hgvsInfo=? WHERE variantAnnotationID=?" , (warning_message,variantAnnotationID))

	# STEP 5 : c_to_p
	try:
		p = am37.c_to_p(c)
	except Exception as e:
	# except hgvs.exceptions.HGVSError as e:
		warning_message = 'warning : %s' % str(e)
		print "\t- %s" % warning_message
		db_cur.execute("UPDATE VariantAnnotation SET hgvs='warning', hgvsInfo=? WHERE variantAnnotationID=?" , (warning_message,variantAnnotationID))
		continue
	p_pos = str(p).split(':')[-1]
	print "\t- protein description : %s" % str(p)
	db_cur.execute("UPDATE VariantAnnotation SET proteinDescription=? WHERE variantAnnotationID=?" , (p_pos,variantAnnotationID))

# MERGE VARIANTMETRICS DUPLICATES (WITH SAME VARIANT) AFTER HGVS CORRECTION
print "- Merge variantMetrics duplicates"
db_cur.execute("SELECT analysis,variant, COUNT(*) c FROM VariantMetrics GROUP BY analysis,variant HAVING c > 1")
db_vm_duplicates = db_cur.fetchall()
for db_vm_dup in db_vm_duplicates:
	db_cur.execute("SELECT * FROM VariantMetrics WHERE analysis='%s' and variant='%s'" % (db_vm_dup['analysis'],db_vm_dup['variant']))
	db_vms_to_merge = db_cur.fetchall()
	vmids = []
	depths = []
	aos = []
	calls = []
	for db_vm_to_merge in db_vms_to_merge:
		vmids.append(db_vm_to_merge['variantMetricsID'])
		depths.append(db_vm_to_merge['positionReadDepth'])
		aos.append(db_vm_to_merge['variantReadDepth'])
		if ' / ' in db_vm_to_merge['call']:
			cc = db_vm_to_merge['call'].split(' / ')
		else:
			cc = db_vm_to_merge['call'].split('/')
		for c in cc :
			calls.append(c)
	depth = int(numpy.round(numpy.mean(depths)))
	ao = int(numpy.round(numpy.mean(aos)))
	calls = list(set(calls))
	call = '/'.join(calls)

	# update de la premiere entree et delete des autres
	db_cur.execute("UPDATE VariantMetrics SET positionReadDepth=%s, variantReadDepth=%s, call='%s' WHERE variantMetricsID='%s'" % (depth,ao,call,vmids[0]))
	for vmid in vmids[1:]:
		db_cur.execute("DELETE FROM VariantMetrics WHERE VariantMetricsID='%s'" % vmid)

db_con.commit()
db_con.close()

print "\t - [%s] done." % (time.strftime("%H:%M:%S"))
