#!/usr/bin/python
from optparse import OptionParser
from datetime import date
import subprocess
import sqlite3
import urllib2
import zipfile
import pysam
import json
import uuid
import time
import sys
import os

# THIS SCRIPT RUN ANNOVAR AND VEP FOR VARIANTS IN VariantBase.db
# USAGE : python annotate_variantbase.step1.py --option

# annovar input format : chr1   115256572   115256572   C   T
# vep input format     :    1   115256572   115256572   C/T   +

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

### GATHERING PARAMETERS ############################################################

parser = OptionParser()
parser.add_option('-n', '--new',	help="annonate new variants only", 	dest='new', 	default=False, action='store_true')
parser.add_option('-f', '--full',	help="annonate full database", 		dest='full', 	default=False, action='store_true')
# parser.add_option('-u', '--update',	help="update old variants", 		dest='update', 	default=False, action='store_true')
parser.add_option('-o', '--output-folder',	help="output folder", 		dest='output')
(options, args) = parser.parse_args()

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))
	
db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

if not os.path.isdir('%s/annovar' % options.output):
	subprocess.call(['mkdir', '%s/annovar' % options.output])
if not os.path.isdir('%s/vep' % options.output):
	subprocess.call(['mkdir', '%s/vep' % options.output])

###################################################################################################################################

# if options.update:
	# exit()
	#db_cur.execute("SELECT * FROM VariantAnnotation WHERE lastUpdate is NULL")  # option pas prete # ici faire test WHERE lastUpdate > 6 mois
if options.new:
	db_cur.execute("""SELECT DISTINCT chromosome,genomicStart,genomicStop,referenceAllele,alternativeAllele,variantType FROM Variant 
	INNER JOIN VariantAnnotation ON VariantAnnotation.variant=Variant.variantID 
	WHERE (lastUpdate is NULL) OR (region is NULL)""")
if options.full:
	db_cur.execute("""SELECT DISTINCT chromosome,genomicStart,genomicStop,referenceAllele,alternativeAllele,variantType FROM Variant 
	INNER JOIN VariantAnnotation ON VariantAnnotation.variant=Variant.variantID """)

db_variants = db_cur.fetchall()
if not db_variants:
	print "- 0 new variants to process"
	exit()

# WRITE INPUT FILES FOR ANNOVAR AND VEP
print "- Get all variants data..."
variant_data = []
for db_variant in db_variants:
	variant_data.append([db_variant['chromosome'],db_variant['genomicStart'],db_variant['genomicStop'],db_variant['referenceAllele'],db_variant['alternativeAllele'],db_variant['variantType']])
sorted_variant_data = sorted(sorted(sorted(variant_data, key = lambda x : int(x[2])), key = lambda x : int(x[1])), key = lambda x : int(x[0].replace('chr','').replace('X','23').replace('Y','24')))
print "\t- %s variants to process..." % len(sorted_variant_data)

print "- Generate Annovar input..."	
with open('%s/annovar/annovar_input.tsv' % options.output,'w') as annovar_input:
	for v in sorted_variant_data:
		line = '%s\t%s\t%s\t%s\t%s\n' % (v[0],v[1],v[2],v[3],v[4])
		if (v[5] == 'INS') or (v[5] == 'DUP'):
			line = '%s\t%s\t%s\t%s\t%s\n' % (v[0],v[1],v[1],v[3],v[4])
		annovar_input.write(line)

print "- Generate VEP input..."	
with open('%s/vep/vep_input.tsv' % options.output,'w') as vep_input:
	for v in sorted_variant_data:
		line = '%s\t%s\t%s\t%s/%s\t+\n' % (v[0].replace('chr',''),v[1],v[2],v[3],v[4])
		if (v[5] == 'INS') or (v[5] == 'DUP'):
			line = '%s\t%s\t%s\t%s/%s\t+\n' % (v[0].replace('chr',''),v[2],v[1],v[3],v[4])
		vep_input.write(line)

print "- [%s] table_annovar.pl ..." % (time.strftime("%H:%M:%S"))
cmd = subprocess.Popen([
	'perl', '%s/variantAnnotation/annovar/table_annovar.pl' % pipeline_folder, 
	'%s/annovar/annovar_input.tsv' % options.output, '%s/variantAnnotation/annovar/humandb/' % pipeline_folder,
	'-buildver', 'hg19',
	'-out', '%s/annovar/annovar' % options.output,
	'-protocol', 'refGeneWithVer,cosmic92,avsnp150,intervar_20180118,clinvar_20200316,nci60,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eur,gnomad211_exome,exac03,dbnsfp35a',
	'-operation', 'g,f,f,f,f,f,f,f,f,f,f,f',
	'-argument', '--hgvs,,,,,,,,,,,',
	'-nastring', '.',
	'-polish',
	'-xref', '%s/variantAnnotation/annovar/example/gene_xref.txt' % pipeline_folder,
	], stdout=open('%s/annovar/table_annovar.stdout.txt' % options.output,'w'), stderr=open('%s/annovar/table_annovar.stderr.txt' % options.output,'w'))
cmd.communicate()

print "- [%s] run_vep ..." % (time.strftime("%H:%M:%S"))
cmd = subprocess.Popen([
	'perl', '%s/variantAnnotation/VEP/ensembl-vep/vep' % pipeline_folder,
	'--offline',
	'--dir_cache', '%s/variantAnnotation/VEP/cache/' % pipeline_folder,
	'--force_overwrite',
	'--refseq',
	'--use_transcript_ref',
	'--numbers',
	'--hgvs',
	'--hgvsg',
	'--no_escape',
	'--check_existing', # retire des "existing variants" les variants connus qui ne correspondent pas aux alleles mutes (par defaut utile seulement les coordonnees)
	# For some data sources (COSMIC, HGMD), Ensembl is not licensed to redistribute allele-specific data, so VEP will report the existence of co-located variants with unknown alleles without carrying out allele matching.
	'--variant_class',
	'--sift', 'p',
	'--polyphen', 'p',
	'--af_1kg',
	'--af',
	'--freq_pop', '1KG_ALL',
	'--freq_pop', 'gnomAD',
	'--af_esp',
	'--symbol',
	'--tab',
	'--pubmed',
	'--fasta', '%s/variantAnnotation/VEP/cache/homo_sapiens/96_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz' % pipeline_folder,
	'--exclude_predicted', # RETIRE LES TRANSCRITS PREDITS "XM" ou "XR"
	'--no_stats',
	'--verbose',
	'--input_file', '%s/vep/vep_input.tsv' % options.output,
	'--format', 'ensembl',
	'--output_file', '%s/vep/vep_output.tsv' % options.output
	], stdout=open('%s/vep/vep.stdout.txt' % options.output,'w'), stderr=open('%s/vep/vep.stderr.txt' % options.output,'w'))
cmd.communicate()

print "[%s] done." % (time.strftime("%H:%M:%S"))
