#!/usr/bin/python
import os
import sys
import csv
import json
import uuid
import time
import zipfile
import sqlite3
from datetime import date
from optparse import OptionParser
from pprint import pprint
from numpy import mean

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
parser.add_option('-v', '--vcf',		help="VCF file (one or more : --vcf 1.vcf --vcf 2.vcf)", action="append", dest='vcfs')
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
genome_build = 'hg19'
vc_tool = '?'

db_cur.execute("SELECT * FROM Gene")
db_genes = db_cur.fetchall()
gene_data = {}
for db_gene in db_genes:
	gene_data[db_gene['geneID']] = {
		'chromosome': db_gene['chromosome'],
		'transcriptionStart': db_gene['transcriptionStart'],
		'transcriptionStop': db_gene['transcriptionStop']
	}

abl1_c2g = {}
if options.abl1 == 'yes':
	abl1_cdna2genomic = global_param['abl1_cdna2genomic']
	abl1_cdna2genomic_file = open(abl1_cdna2genomic,'r')
	abl1_cdna2genomic_reader = csv.reader(abl1_cdna2genomic_file,delimiter='\t')
	abl1_cdna2genomic_reader.next() # header
	for line in abl1_cdna2genomic_reader:
		abl1_c2g[int(line[2])] = int(line[3])

mutect2_nocall = ['base_qual','low_allele_frac','map_qual','slippage','strand_bias','strict_strand','weak_evidence']

#   __        __   __          __       __   ___  __            ___  __  
#  |__)  /\  |__) /__` | |\ | / _`     |__) |__  /__` |  | |     |  /__` 
#  |    /~~\ |  \ .__/ | | \| \__>     |  \ |___ .__/ \__/ |___  |  .__/ 

variants = {}
mutect2_count = 0
lofreq_count = 0
varscan_count = 0
vardict_count = 0
# deepvariant_count = 0

for vcf_path in options.vcfs:
	vcf = open(vcf_path,'r')
	for line in vcf:
		if line.startswith('#'):
			continue
		line = line.replace('\n','').split('\t')
		chromosome = line[0]
		position = int(line[1])
		oref = line[3]
		alts = line[4].split(',')
		if 'mutect' in vcf_path:
			vc_name = 'mutect2'
			filters = line[6].split(';')
			pass_variant = False
			for no_call in mutect2_nocall:
				if no_call in filters:
					pass_variant = True
			if pass_variant:
				continue
			pos_cov = int(line[9].split(':')[3])
			var_covs = line[9].split(':')[1].split(',')[1:] # premiere valeur est ref count
			var_covs = [int(v) for v in var_covs]
			mutect2_count += 1
		elif 'lofreq' in vcf_path:# 1 seul variant par ligne avec lofreq
			vc_name = 'lofreq'
			dp4 = line[7].split('DP4=')[-1].split(';')[0].split(',') 
			pos_cov = int(dp4[0])+int(dp4[1])+int(dp4[2])+int(dp4[3])
			var_covs = [int(dp4[2])+int(dp4[3])]
			lofreq_count += 1
		elif 'varscan' in vcf_path: #VarScan calls the most-observed variant
			vc_name = 'varscan2'
			filters = line[6]
			if filters != 'PASS':
				continue
			pos_cov = int(line[9].split(':')[2])
			var_covs = [int(line[9].split(':')[5])]
			varscan_count += 1
		elif 'vardict' in vcf_path:
			vc_name = 'vardict'
			filters = line[6]
			if filters != 'PASS':
				continue
			pos_cov = int(line[9].split(':')[1])
			var_covs = line[9].split(':')[3].split(',')[1:] # premiere valeur est ref count
			var_covs = [int(v) for v in var_covs]
			vardict_count += 1
		# elif 'deepvariant' in vcf_path:
			# vc_name = 'deepvariant'
			# filters = line[6]
			# if filters != 'PASS':
				# continue
			# pos_cov = int(line[9].split(':')[2])
			# var_covs = line[9].split(':')[3].split(',')[1:] # premiere valeur est ref count
			# var_covs = [int(v) for v in var_covs]
			# deepvariant_count += 1

		for a in range(len(alts)): # si multiallelic
			ref = oref
			alt = alts[a]
			var_cov = var_covs[a]
			# left-align
			while ref[0] == alt[0]:
				ref = ref[1:]
				alt = alt[1:]
				position += 1
				if ref == '':
					ref = '-'
				if alt == '':
					alt = '-'
			# right-align
			while ref[-1] == alt[-1]:
				ref = ref[:-1]
				alt = alt[:-1]
				if ref == '':
					ref = '-'
				if alt == '':
					alt = '-'
			#Genomic Description (do not forget position has been modified with left-alignment)
			if ref == '-': # INS
				variant_type = 'INS'
				start = position - 1
				stop = position
				genomicDescription = '%s:g.%s_%sins%s' % (chromosome,start,stop,alt)
			elif alt == '-': # DEL
				variant_type = 'DEL'
				if len(ref) > 1:
					start = position
					stop = position + (len(ref)-1)
					genomicDescription = '%s:g.%s_%sdel%s' % (chromosome,start,stop,ref)
				else:
					start = position
					stop = position
					genomicDescription = '%s:g.%sdel%s' % (chromosome,start,ref)
			elif len(ref) > 1 or len(alt) > 1: # DELINS
				variant_type = 'DELINS'
				if len(ref) > 1:
					start = position
					stop = position + (len(ref)-1)
					genomicDescription = '%s:g.%s_%sdelins%s' % (chromosome,start,stop,alt)
				else:
					start = position
					stop = position
					genomicDescription = '%s:g.%sdelins%s' % (chromosome,start,alt)
			else: # SNP
				variant_type = 'SNP'
				start = position
				stop = position
				genomicDescription = '%s:g.%s%s>%s' % (chromosome,start,ref,alt)
				
			gene = ''
			for g in gene_data.keys():
				if chromosome == gene_data[g]['chromosome'] and ((gene_data[g]['transcriptionStart']-5000)<start<(gene_data[g]['transcriptionStop']+5000)):
					gene = str(g)
			
			variant = '%s:%s-%s:%s>%s' % (chromosome,start,stop,ref,alt)
			if variant not in variants.keys():
				variants[variant] = {
					'chromosome':chromosome,
					'start':start,
					'stop':stop,
					'ref':ref,
					'alt':alt,
					'variant_type':variant_type,
					'gene':gene,
					'pos_cov':[pos_cov],
					'var_cov':[var_cov],
					'call':[vc_name],
					'genomicDescription':genomicDescription
				}
			elif vc_name not in variants[variant]['call']: # evite lofreq lignes en double...
				variants[variant]['pos_cov'].append(pos_cov)
				variants[variant]['var_cov'].append(var_cov)
				variants[variant]['call'].append(vc_name)

print " variants called by Mutect2     : %s" % mutect2_count
print " variants called by Lofreq      : %s" % lofreq_count
print " variants called by VarScan     : %s" % varscan_count
print " variants called by VarDict     : %s" % vardict_count
# print " variants called by Deepvariant : %s" % deepvariant_count
print " - TOTAL : %s variants" % len(variants)
x = 0
y = 0
z = 0
for variant in variants:
	if len(variants[variant]['call']) >= 2:
		x+=1
	if len(variants[variant]['call']) >= 3:
		y+=1
	if len(variants[variant]['call']) == 4:
		z+=1
print " - found by 4 callers  : %s variants" % z
print " - found by 3+ callers : %s variants" % y
print " - found by 2+ callers : %s variants" % x

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
	if (len(variants[variant]['call']) >= 2) or (variants[variant]['gene'] == 'CEBPA' and variants[variant]['variant_type'] != 'SNP'):
		variant_id = variant
		pos_cov = int(mean(variants[variant]['pos_cov']))
		var_cov = int(mean(variants[variant]['var_cov']))
		call = '/'.join(variants[variant]['call'])
		db_cur.execute("SELECT * FROM Variant WHERE variantID='%s'" % variant_id)
		db_variant = db_cur.fetchone()
		if db_variant is None:
			try:
				#print "- [Variant] : adding %s in DB" % variant_id
				db_cur.execute("INSERT INTO Variant (variantID, genomeBuild, chromosome, genomicStart, genomicStop, referenceAllele, alternativeAllele, variantType, gene, genomicDescription) VALUES ('%s','%s','%s',%s, %s,'%s','%s','%s','%s','%s')" % (variant_id, genome_build, variants[variant]['chromosome'], variants[variant]['start'], variants[variant]['stop'], variants[variant]['ref'], variants[variant]['alt'],variants[variant]['variant_type'],variants[variant]['gene'],variants[variant]['genomicDescription']))
				new_var_count += 1
			except Exception as e:
				print "\t - warning (VARIANT table)** %s"%e
		elif db_variant['hgvs'] == 'no':
			variant_id = db_variant['hgvsInfo']

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
				db_cur.execute("INSERT INTO VariantMetrics (variantMetricsID, variant, analysis, positionReadDepth, variantReadDepth, variantCallingTool, call) VALUES ('%s', '%s', '%s', %s, %s, '%s', '%s')" % (variantmetrics_id, variant_id, analysis_id, pos_cov, var_cov, vc_tool, call))
				new_vm_count += 1
			except Exception as e:
				print "\t - warning (VARIANTMETRICS table)** %s"%e

print " - [%s] %s variants (%s new entries, %s occurences added)" % (time.strftime("%H:%M:%S"),z,new_var_count,new_vm_count)

db_con.commit()
db_con.close()
