#!/usr/bin/python
import os
import csv
import json
import uuid
import time
import zipfile
import sqlite3
from datetime import date
from optparse import OptionParser

# Using variant results file alleles.txt, this script update VariantBase insertings new variants and variantMetrics
# USAGE : python insert_db_variants.py --analysis analysisID --variants /.../alleles.xls

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

### GATHERING PARAMETERS ############################################################

parser = OptionParser()
parser.add_option('-a', '--analysis',	help="DB analysisID", dest='analysis')
parser.add_option('-v', '--variants',	help="TVC alleles.xls", dest='variants')
parser.add_option('-z', '--abl1',		help="abl1 conversion", dest='abl1', default='no')
(options, args) = parser.parse_args()

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

analysis_id = options.analysis
vc_tool = 'TVC'

vc_xls_path = options.variants
if os.path.isfile(vc_xls_path):
	vc_xls_file = open(vc_xls_path,'r')
else:
	intermediate_folder = os.path.dirname(os.path.dirname(vc_xls_path))
	if os.path.isfile('%s.zip' % intermediate_folder):
		archive = zipfile.ZipFile('%s.zip' % intermediate_folder, 'r')
		if 'tvc_de_novo/alleles.xls' in archive.namelist():
			vc_xls_file = archive.open('tvc_de_novo/alleles.xls')	

vc_xls_reader = csv.DictReader(vc_xls_file,delimiter='\t')

vc_hotspot_xls_reader = False
if os.path.isfile(vc_xls_path.replace('tvc_de_novo','tvc_only_hotspot')):
	vc_hotspot_xls_file = open(vc_xls_path.replace('tvc_de_novo','tvc_only_hotspot'),'r')
	vc_hotspot_xls_reader = csv.DictReader(vc_hotspot_xls_file,delimiter='\t')
else:
	intermediate_folder = os.path.dirname(os.path.dirname(vc_xls_path))
	if os.path.isfile('%s.zip' % intermediate_folder):
		archive = zipfile.ZipFile('%s.zip' % intermediate_folder, 'r')
		if 'tvc_only_hotspot/alleles.xls' in archive.namelist():
			vc_xls_file = archive.open('tvc_only_hotspot/alleles.xls')
			vc_hotspot_xls_reader = csv.DictReader(vc_hotspot_xls_file,delimiter='\t')

abl1_c2g = {}
# /DATA/work/variantAnnotation/ABL1_NM_005157_cdna_genomic.tsv
# http://grch37.ensembl.org/biomart
# Filters : RefSeq mRNA IDs : NM_005157
# Attributes : Chromosome/scaffold name, Exon rank in transcript, Genomic coding start, Genomic coding end, CDS start, CDS end
if options.abl1 == 'yes':
	abl1_cdna2genomic = global_param['abl1_cdna2genomic']
	abl1_cdna2genomic_file = open(abl1_cdna2genomic,'r')
	abl1_cdna2genomic_reader = csv.reader(abl1_cdna2genomic_file,delimiter='\t')
	abl1_cdna2genomic_reader.next() # header
	for line in abl1_cdna2genomic_reader:
		abl1_c2g[int(line[2])] = int(line[3])

#   __        __   __          __     ___       __      __   ___  __            ___  __  
#  |__)  /\  |__) /__` | |\ | / _`     |  \  / /  `    |__) |__  /__` |  | |     |  /__` 
#  |    /~~\ |  \ .__/ | | \| \__>     |   \/  \__,    |  \ |___ .__/ \__/ |___  |  .__/ 

variants = {}

for line in vc_xls_reader:
	call = line['Allele Call']
	if call != 'Heterozygous' and call != 'Homozygous':
		continue
	chromosome = line['Chrom']
	position = int(line['Position'])
	if options.abl1 == 'yes':
		chromosome = 'chr9'
		position = abl1_c2g[position-192]
	ref = line['Ref']
	alt = line['Variant']
	variant_type = line['Type']
	freq = line['Frequency']
	pos_cov = int(line['Coverage'])
	var_cov = int(line['Allele Cov'])

	if variant_type == 'SNP':
		start = position
		stop = position
	elif variant_type == 'DEL':
		start = position
		stop = position+len(ref)-1
	elif variant_type == 'INS':
		start = position-1
		stop = position
	elif variant_type == 'COMPLEX' or variant_type == 'MNP':
		start = position
		stop = position+len(ref)-1

	variant = '%s:%s-%s:%s>%s' % (chromosome,start,stop,ref,alt)
	if variant not in variants.keys():
		variants[variant] = {
			'chromosome':chromosome,
			'start':start,
			'stop':stop,
			'ref':ref,
			'alt':alt,
			'variant_type':variant_type,
			'freq':freq,
			'pos_cov':pos_cov,
			'var_cov':var_cov,
			'call':'tvc_de_novo'
		}

if vc_hotspot_xls_reader:
	for line in vc_hotspot_xls_reader:
		call = line['Allele Call']
		if call != 'Heterozygous' and call != 'Homozygous':
			continue
		chromosome = line['Chrom']
		position = int(line['Position'])
		if options.abl1 == 'yes':
			chromosome = 'chr9'
			position = abl1_c2g[position-192]
		ref = line['Ref']
		alt = line['Variant']
		variant_type = line['Type']
		freq = line['Frequency']
		pos_cov = int(line['Coverage'])
		var_cov = int(line['Allele Cov'])

		if variant_type == 'SNP':
			start = position
			stop = position
		elif variant_type == 'DEL':
			start = position
			stop = position+len(ref)-1
		elif variant_type == 'INS':
			start = position-1
			stop = position
		elif variant_type == 'COMPLEX' or variant_type == 'MNP':
			start = position
			stop = position+len(ref)-1

		variant = '%s:%s-%s:%s>%s' % (chromosome,start,stop,ref,alt)
		if variant not in variants.keys():
			variants[variant] = {
				'chromosome':chromosome,
				'start':start,
				'stop':stop,
				'ref':ref,
				'alt':alt,
				'variant_type':variant_type,
				'freq':freq,
				'pos_cov':pos_cov,
				'var_cov':var_cov,
				'call':'tvc_hotspot'
			}
		else:
			variants[variant]['call'] = 'tvc_de_novo/tvc_hotspot'

# IF ANALYSIS PREVIOUSLY DONE, DELETE OLD VARIANTMETRICS FIRST
db_cur.execute("SELECT variantMetricsID FROM VariantMetrics WHERE analysis='%s'" % analysis_id)
db_vms = db_cur.fetchall()
if db_vms:
	print " - [%s] analysisID already in DB : cleaning previous variantmetrics" % (time.strftime("%H:%M:%S"))
	for db_vm in db_vms:
		db_cur.execute("DELETE FROM VariantMetrics WHERE variantMetricsID='%s'" % db_vm['variantMetricsID'])

#   ___                     __                __              ___    ___       __        ___ 
#  |__  | |    |    | |\ | / _`    \  /  /\  |__) |  /\  |\ |  |      |   /\  |__) |    |__  
#  |    | |___ |___ | | \| \__>     \/  /~~\ |  \ | /~~\ | \|  |      |  /~~\ |__) |___ |___ 

new_var_count = 0
new_vm_count = 0

for variant in variants:
	variant_id = variant

	# CREATE NEW ENTRY IF VARIANT DOESNT EXISTS
	db_cur.execute("SELECT * FROM Variant WHERE variantID='%s'" % variant_id)
	db_variant = db_cur.fetchone()
	if db_variant is None:
		try:
			print "- [Variant] : new entry %s in DB" % variant_id
			db_cur.execute("INSERT INTO Variant (variantID, chromosome, genomicStart, genomicStop, referenceAllele, alternativeAllele) VALUES ('%s','%s',%s, %s,'%s','%s')" % (variant_id, variants[variant]['chromosome'], variants[variant]['start'], variants[variant]['stop'], variants[variant]['ref'], variants[variant]['alt']))
			new_var_count += 1
		except Exception as e:
			print "\t - warning (VARIANT table)** %s"%e

	#   ___                     __                __              ___        ___ ___  __     __   __     ___       __        ___ 
	#  |__  | |    |    | |\ | / _`    \  /  /\  |__) |  /\  |\ |  |   |\/| |__   |  |__) | /  ` /__`     |   /\  |__) |    |__  
	#  |    | |___ |___ | | \| \__>     \/  /~~\ |  \ | /~~\ | \|  |   |  | |___  |  |  \ | \__, .__/     |  /~~\ |__) |___ |___ 

	db_cur.execute("SELECT variantMetricsID FROM VariantMetrics WHERE analysis='%s' and variant='%s'" % (analysis_id,variant_id))
	if db_cur.fetchone() is None: 
		exists = True
		while exists != None:
			random_uuid = uuid.uuid1()
			variantmetrics_id = 'M-'+random_uuid.hex[:8]
			db_cur.execute("SELECT variantMetricsID FROM VariantMetrics WHERE variantMetricsID='%s'" % (variantmetrics_id))
			exists = db_cur.fetchone()
		try:
			#print "- [VariantMetrics] : adding %s in DB (%s <-> %s)" % (variantmetrics_id,variant_id,analysis_id)
			db_cur.execute("INSERT INTO VariantMetrics (variantMetricsID, variant, analysis, positionReadDepth, variantReadDepth, variantCallingTool, call) VALUES ('%s', '%s', '%s', %s, %s, '%s', '%s')" % (variantmetrics_id, variant_id, analysis_id, variants[variant]['pos_cov'], variants[variant]['var_cov'], vc_tool, variants[variant]['call']))
			new_vm_count += 1
		except Exception as e:
			print "\t - warning (VARIANTMETRICS table)** %s"%e

print " - [%s] %s variants (%s new Variant entries, %s occurences added)" % (time.strftime("%H:%M:%S"),len(variants.keys()),new_var_count,new_vm_count)

db_con.commit()
db_con.close()
