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
# import docker
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
parser.add_option('-n', '--new',	help="annonate new variants only", 	dest='new', 	default=False, action='store_true')
parser.add_option('-f', '--full',	help="annonate full database", 		dest='full', 	default=False, action='store_true')
parser.add_option('-u', '--update',	help="update old variants", 		dest='update', 	default=False, action='store_true')
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

# try:
	# hdp = hgvs.dataproviders.uta.connect()
# except:
	# client = docker.from_env()
	# container = client.containers.get(global_param['uta_container_id'])
	# print "- uta container not responding : %s" % container
	# print "- trying to restart docker uta container..."
	# container.restart()
	# time.sleep(10)
	# hdp = hgvs.dataproviders.uta.connect()
hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
hn = hgvs.normalizer.Normalizer(hdp)
hv = hgvs.validator.Validator(hdp)
am37 = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37')
am37no = hgvs.assemblymapper.AssemblyMapper(hdp, normalize=False, assembly_name='GRCh37') # for c to g conversion when c already hgvs

###################################################################################################################################

# FULL HGVS c. / p. RE-ANNOTATION

if options.update:
	exit()
	# faire WHERE hgvs is NULL OR lastUpdate > 6 mois ?
	# db_cur.execute("SELECT * FROM Variant WHERE lastUpdate > date('now', 'start of day','-5 days')")
	# UTILISER DATETIME https://stackoverflow.com/questions/1975737/sqlite-datetime-comparison
if options.new:
	# db_cur.execute("SELECT * FROM Variant WHERE hgvs is NULL")
	db_cur.execute("SELECT * FROM Variant WHERE (hgvs is NULL) OR (hgvs='yes' AND region is NULL)")
if options.full:
	db_cur.execute("SELECT * FROM Variant")
db_variants = db_cur.fetchall()
print "- %s new variants to process" % len(db_variants)

lastUpdate = '%s' % (time.strftime('%d/%m/%Y'))
z=0
for db_variant in db_variants:
	z+=1
	variantID = db_variant['variantID']
	genomeBuild = db_variant['genomeBuild']
	chromosome = db_variant['chromosome']
	genomicStart = db_variant['genomicStart']
	genomicStop = db_variant['genomicStop']
	ref = db_variant['referenceAllele']
	alt = db_variant['alternativeAllele']
	gene = db_variant['gene']
	variantType = db_variant['variantType']
	genomicDescription = db_variant['genomicDescription'].split(':')[-1]
	db_cur.execute("SELECT * FROM Gene WHERE geneID='%s'" % gene)
	db_gene = db_cur.fetchone()
	if not db_gene:
		print "- %s: WARNING : Gene not found in DB (%s:%s-%s:%s>%s)" % (z,chromosome,genomicStart,genomicStop,ref,alt)
	transcript = "%s.%s" % (db_gene['transcript'],db_gene['transcriptVersion'])
	nc = db_gene['NC']
	strand = db_gene['strand']

	# STEP 1 : PARSE VARIANT (g)
	g0 = hp.parse_hgvs_variant('%s:%s' % (nc,genomicDescription))
	print "- %s: %s (%s:%s-%s:%s>%s)" % (z,g0,chromosome,genomicStart,genomicStop,ref,alt)
	
	# STEP 2 : VALIDATE (g)
	try:
		hv.validate(g0)
	except hgvs.exceptions.HGVSError as e:
		error_message = 'error : %s' % str(e)
		print "\t- %s" % error_message
		db_cur.execute("UPDATE Variant SET hgvs='error', hgvsInfo='%s' WHERE variantID='%s'" % (error_message,variantID))
		continue
		
	# STEP 3 : g to c (with normalization)
	try:
		c = am37.g_to_c(g0,transcript)
	except hgvs.exceptions.HGVSError as e:
		error_message = 'error : %s' % str(e)
		print "\t- %s" % error_message
		db_cur.execute("UPDATE Variant SET hgvs='error', hgvsInfo='%s' WHERE variantID='%s'" % (error_message,variantID))
		continue
		
	g = am37no.c_to_g(c) # Get coherent g (if not already) from normalyzed c (3' rule)
	# why not normalize g first? normalize(g) is useless, it just right-shift in the genomic sens (no 3' rule)
	
	# STEP 4 : is the variant HGVS? if not create HGVS entry in db
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
			#pas besoin de nouvelle entree, mais update genomic description et le reste
			genomicDescription = '%s:%s' % (chromosome,str(g).split(':')[-1])
			db_cur.execute("UPDATE Variant SET hgvs='yes', genomicDescription='%s', variantType='%s', lastUpdate='%s' WHERE variantID='%s'" % (genomicDescription,variantType,lastUpdate,variantID))
		else:
			bad_variantID = variantID
			variantID = '%s:%s-%s:%s>%s' % (chromosome,start,end,ref,alt)
			genomicDescription = '%s:%s' % (chromosome,str(g).split(':')[-1])
			db_cur.execute("SELECT * FROM Variant WHERE variantID='%s'"%variantID)
			if db_cur.fetchone() is None:
				print "\t\t- new entry created in db for %s" % variantID
				db_cur.execute("INSERT INTO Variant (variantID, genomeBuild, chromosome, genomicStart, genomicStop, referenceAllele, alternativeAllele, variantType, gene, genomicDescription) VALUES ('%s','%s','%s',%s, %s,'%s','%s','%s','%s','%s')" % (variantID,genomeBuild,chromosome,start,end,ref,alt,variantType,gene,genomicDescription))
			db_cur.execute("UPDATE Variant SET hgvs='no', hgvsInfo='%s', lastUpdate='%s' WHERE variantID = '%s'" % (variantID,lastUpdate,bad_variantID))
			db_cur.execute("SELECT * FROM VariantMetrics WHERE variant = '%s'" % bad_variantID)
			varmetrics = db_cur.fetchall()
			for varmetric in varmetrics:
				VariantMetricsID = varmetric['variantMetricsID']
				db_cur.execute("UPDATE VariantMetrics SET variant = '%s' WHERE VariantMetricsID = '%s'" % (variantID,VariantMetricsID))
			
	# STEP 5 : VALIDATE(c)
	c_pos = str(c).split(':')[-1]
	print "\t- %s" % c_pos
	db_cur.execute("UPDATE Variant SET transcriptDescription='%s', hgvs='yes' WHERE variantID='%s'" % (c_pos,variantID))
	
	try:
		hv.validate(c)
	except hgvs.exceptions.HGVSError as e:
		warning_message = 'warning : %s' % str(e)
		print "\t- %s" % warning_message
		db_cur.execute("UPDATE Variant SET hgvs='warning', hgvsInfo='%s' WHERE variantID='%s'" % (warning_message,variantID))

	# STEP 5 : c_to_p
	try:
		p = am37.c_to_p(c)
	except hgvs.exceptions.HGVSError as e:
		warning_message = 'warning : %s' % str(e)
		print "\t- %s" % warning_message
		db_cur.execute("UPDATE Variant SET hgvs='warning', hgvsInfo='%s' WHERE variantID='%s'" % (warning_message,variantID))
		continue
		
	p_pos = str(p).split(':')[-1]
	print "\t- %s" % p_pos
	db_cur.execute("UPDATE Variant SET proteinDescription='%s' WHERE variantID='%s'" % (p_pos,variantID))

	
db_con.commit()
db_con.close()

print "\t - [%s] done." % (time.strftime("%H:%M:%S"))
