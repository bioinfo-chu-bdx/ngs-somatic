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
import json
import time

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
	#db_cur.execute("SELECT * FROM Variant WHERE lastUpdate > date('now', 'start of day','-5 days')")
	# UTILISER DATETIME https://stackoverflow.com/questions/1975737/sqlite-datetime-comparison
if options.new:
	db_cur.execute("SELECT * FROM Variant WHERE hgvs is NULL")
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

	# STEP 1 : PARSE
	g0 = hp.parse_hgvs_variant('%s:%s' % (nc,genomicDescription))
	
	print "- %s: %s (%s:%s-%s:%s>%s)" % (z,g0,chromosome,genomicStart,genomicStop,ref,alt)
	# STEP 2 : VALIDATE
	try:
		hv.validate(g0)
	except hgvs.exceptions.HGVSError as e:
		print '\t-G-ERROR : %s' % (e)
		continue
		
	# STEP 3 : g -> c (+ normalize)
	try:
		c = am37.g_to_c(g0,transcript)
	except hgvs.exceptions.HGVSError as e:
		print '\t-G2C-ERROR : %s' % (e)
		continue
		
	try:
		hv.validate(c)
	except hgvs.exceptions.HGVSError as e:
		print '\t-C-WARNING : %s' % (e)
	print '\t-'+str(c).split(':')[-1]
	
	# STEP 5 : c_to_p
	try:
		p = am37.c_to_p(c)
	except hgvs.exceptions.HGVSError as e:
		print '\t-P-ERROR : %s' % (e)
		
	#try:
		#hv.validate(p)
	#except hgvs.exceptions.HGVSError as e:
		#print '\t-P-WARNING : %s' % (e)
	print '\t-'+str(p).split(':')[-1]
		
	g = am37no.c_to_g(c) # get g from good c without norm. 
	
print "\t - [%s] done." % (time.strftime("%H:%M:%S"))
