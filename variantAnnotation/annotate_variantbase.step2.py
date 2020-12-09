#!/usr/bin/python
import os
import sys
import csv
import json
import time
import xlrd
import sqlite3
from optparse import OptionParser

# THIS SCRIPT USE ANNOVAR AND VEP RESULTS FOR REANNOTATE ALL VARIANTS IN VariantBase.db
# USAGE : python reannotate_all_variantbase.step2.py --annovar-results /path/to/annovar_results --vep-results /path/to/vep_results

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

def representsInt(s):
	try: 
		s = int(s)
		return s
	except:
		return s

def representsFloat(s):
	try: 
		s = float(s)
		return s
	except:
		return s

def compare_consequence(region_consequence_list):
	most_severe_region = None
	most_severe_consequence = None
	regions = []
	consequences = []
	for tup in region_consequence_list:
		regions.append(tup[0])
		consequences.append(tup[1])
	## REGION ##
	if 'exonic;splicing' in regions or ('exonic' in regions and 'splicing' in regions):
		most_severe_region = 'exonic;splicing'
	elif 'splicing' in regions:
		most_severe_region = 'splicing'
	elif 'exonic' in regions:
		most_severe_region = 'exonic'
	elif 'ncRNA' in regions:
		most_severe_region = 'ncRNA'
	elif 'UTR5' in regions:
		most_severe_region = 'UTR5'
	elif 'UTR3' in regions:
		most_severe_region = 'UTR3'
	elif 'intronic' in regions:
		most_severe_region = 'intronic'
	elif 'upstream' in regions:
		most_severe_region = 'upstream'
	elif 'downstream' in regions:
		most_severe_region = 'downstream'
	## consequence ##
	if 'frameshift' in consequences:
		most_severe_consequence = 'frameshift'
	elif 'stopgain' in consequences:
		most_severe_consequence = 'stopgain'
	elif 'stoploss' in consequences:
		most_severe_consequence = 'stoploss'
	elif 'startloss' in consequences:
		most_severe_consequence = 'startloss'
	elif 'splicing' in consequences:
		most_severe_consequence = 'splicing'
	elif 'nonframeshift' in consequences:
		most_severe_consequence = 'nonframeshift'
	elif 'missense' in consequences:
		most_severe_consequence = 'missense'
	elif 'synonymous' in consequences:
		most_severe_consequence = 'synonymous'
	elif '.' in consequences:
		most_severe_consequence = '.'
	return(most_severe_region,most_severe_consequence)

##########################################################################################################################################
parser = OptionParser()
parser.add_option('-a', '--annovar-results', 		help="Annovar results file (..._multianno.txt)", 	dest='annovar_results')
parser.add_option('-v', '--vep-results', 			help="VEP results file (optionnal)", 				dest='vep_results', default=False)
(options, args) = parser.parse_args()

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

with open('%s/variantAnnotation/vep2annovar.json' % pipeline_folder, 'r') as g:
	vep2annovar = json.load(g)

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

if not os.path.isfile(options.annovar_results):
	print "- 0 new variant annotations to process"
	exit()

annovar_file = open(options.annovar_results,'r')
annovar_reader = csv.DictReader(annovar_file,delimiter='\t')

vep_file = open(options.vep_results,'r')
vep_reader = csv.reader(vep_file,delimiter='\t')

# custom_false_positives = global_param['False_positives']
# custom_drug_sensitivity = global_param['SBT_sensitivity']
# custom_highlight = global_param['LAM_hotspot_highlight']
# custom_tp53_umd = global_param['TP53_UMD_variants_EU']
# custom_lymphome_controls = global_param['LymphomeT_controls']

annotation = {}

############################
# PARSE ANNOTATION RESULTS #
############################
#                   __             __  
#    /\  |\ | |\ | /  \ \  /  /\  |__) 
#   /~~\ | \| | \| \__/  \/  /~~\ |  \ 

print "- [%s] parsing annovar results ..." % (time.strftime("%H:%M:%S"))
z=0
rows=len(list(annovar_reader))
annovar_file.seek(0)
for line in annovar_reader:
	z+=1
	print "\r\t - %s/%s" % (z,rows),
	variantID = '%s:%s-%s:%s>%s' % (line['Chr'],line['Start'],line['End'],line['Ref'],line['Alt'])
	if line['Ref'] == '-':
		variantID = '%s:%s-%s:%s>%s' % (line['Chr'],line['Start'],str(int(line['Start'])+1),line['Ref'],line['Alt'])

	db_cur.execute("""SELECT * FROM VariantAnnotation 
	INNER JOIN Variant ON Variant.variantID=VariantAnnotation.variant
	INNER JOIN Transcript ON Transcript.TranscriptID = VariantAnnotation.transcript
	WHERE variant='%s'""" % variantID)
	db_variant_annotations = db_cur.fetchall() # il peut y avoir plusieurs entrees si plusieurs transcrits pour un meme gene
	for db_variant_annotation in db_variant_annotations:
		variantAnnotationID = db_variant_annotation['variantAnnotationID']
		annotation[variantAnnotationID] = {}
		annotation[variantAnnotationID]['transcript'] = db_variant_annotation['transcript'] # FROM insert_db_variants.py
		annotation[variantAnnotationID]['transcript_without_version'] = db_variant_annotation['transcript'].split('.')[0]
		annotation[variantAnnotationID]['gene'] = db_variant_annotation['gene'].split('.')[0]
		annotation[variantAnnotationID]['genomicStart'] = db_variant_annotation['genomicStart'] # FROM VCF / STEP 0
		annotation[variantAnnotationID]['transcriptDescription'] = db_variant_annotation['transcriptDescription'] # FROM STEP 0
		annotation[variantAnnotationID]['proteinDescription'] = db_variant_annotation['proteinDescription'] # FROM STEP 0

		annotation[variantAnnotationID]['commentaire'] = []
		annotation[variantAnnotationID]['userComment'] = []
		annotation[variantAnnotationID]['annoWarning'] = []
		annotation[variantAnnotationID]['exon'] = None
		annotation[variantAnnotationID]['vep_exon'] = None
		annotation[variantAnnotationID]['vep_intron'] = None
		annotation[variantAnnotationID]['region'] = None
		annotation[variantAnnotationID]['consequence'] = None
		annotation[variantAnnotationID]['vep_consequence'] = None
		annotation[variantAnnotationID]['vep_impact'] = None
		annotation[variantAnnotationID]['annovar_transcriptDescription'] = None
		annotation[variantAnnotationID]['vep_transcriptDescription'] = None
		annotation[variantAnnotationID]['proteinDescription'] = None
		annotation[variantAnnotationID]['annovar_proteinDescription'] = None
		annotation[variantAnnotationID]['actionability'] = None
		annotation[variantAnnotationID]['pubmed'] = None
		annotation[variantAnnotationID]['cosmic'] = None
		annotation[variantAnnotationID]['vep_cosmic'] = None
		annotation[variantAnnotationID]['dbsnp'] = None
		annotation[variantAnnotationID]['vep_dbsnp'] = None
		annotation[variantAnnotationID]['intervar'] = None
		annotation[variantAnnotationID]['clinsig'] = None
		annotation[variantAnnotationID]['vep_clinsig'] = None
		annotation[variantAnnotationID]['nci60'] = None
		annotation[variantAnnotationID]['esp'] = None
		annotation[variantAnnotationID]['vep_esp'] = None
		annotation[variantAnnotationID]['1000G_all'] = None
		annotation[variantAnnotationID]['vep_1000G_all'] = None
		annotation[variantAnnotationID]['1000G_eur'] = None
		annotation[variantAnnotationID]['vep_1000G_eur'] = None
		annotation[variantAnnotationID]['gnomad'] = None
		annotation[variantAnnotationID]['exac'] = None
		annotation[variantAnnotationID]['sift'] = None
		annotation[variantAnnotationID]['vep_sift'] = None
		annotation[variantAnnotationID]['polyphen_hvar'] = None
		annotation[variantAnnotationID]['vep_polyphen_hvar'] = None
		annotation[variantAnnotationID]['provean'] = None
		annotation[variantAnnotationID]['highlight'] = 0
		annotation[variantAnnotationID]['pathoUMD'] = None
		annotation[variantAnnotationID]['commentUMD'] = None
		annotation[variantAnnotationID]['vep_diff'] = None

		# PARSE FUNCTIONNAL ANNOTATION
		transcript_found = False
		if not line['AAChange.refGeneWithVer'] == '.':
			aachange = line['AAChange.refGeneWithVer'].split(',') # format GENE:NM:exon:c:p,
			for item in aachange:
				aadesc = item.split(':')
				if annotation[variantAnnotationID]['transcript_without_version'] == aadesc[1].split('.')[0]:
					transcript_found = True
					if annotation[variantAnnotationID]['transcript'] != aadesc[1]: # not exact transcript version?
						annotation[variantAnnotationID]['annoWarning'].append("Annovar : transcript used is %s (instead of %s)" % (aadesc[1],annotation[variantAnnotationID]['transcript']))
					annotation[variantAnnotationID]['region'] = line['Func.refGeneWithVer']
					annotation[variantAnnotationID]['consequence'] = line['ExonicFunc.refGeneWithVer']
					annotation[variantAnnotationID]['exon'] = representsInt(aadesc[2].replace('exon',''))
					annotation[variantAnnotationID]['annovar_transcriptDescription'] = aadesc[3]
					annotation[variantAnnotationID]['annovar_proteinDescription'] = aadesc[4]
					break
		# non-exonic avec c. (cas particulier : aachange est vide et l'equivalent se trouve dans dans genedetail)
		elif line['GeneDetail.refGeneWithVer'] != '.':
			if 'splicing' in line['Func.refGeneWithVer']:
				cchange = line['GeneDetail.refGeneWithVer'].split(';') # format NM:exon:c;
				for item in cchange:
					cdesc = item.split(':')
					if annotation[variantAnnotationID]['transcript_without_version'] == cdesc[0].split('.')[0]:
						transcript_found = True
						if annotation[variantAnnotationID]['transcript'] != cdesc[0]: # not exact transcript version?
							annotation[variantAnnotationID]['annoWarning'].append("Annovar : transcript used is %s (instead of %s)" % (cdesc[0],annotation[variantAnnotationID]['transcript']))
						annotation[variantAnnotationID]['region'] = line['Func.refGeneWithVer']
						annotation[variantAnnotationID]['consequence'] = line['ExonicFunc.refGeneWithVer']
						annotation[variantAnnotationID]['exon'] = representsInt(cdesc[1].replace('exon',''))
						annotation[variantAnnotationID]['annovar_transcriptDescription'] = cdesc[2]
						break
		else : # ncRNA, UTR3, UTR5, upstream, downstream
			cchange = line['GeneDetail.refGeneWithVer'].split(';')  # format NM:c;
			for item in cchange:
				cdesc = item.split(':') 
				if annotation[variantAnnotationID]['transcript_without_version'] == cdesc[0].split('.')[0]:
					transcript_found = True
					if annotation[variantAnnotationID]['transcript'] != cdesc[0]: # not exact transcript version?
						annotation[variantAnnotationID]['annoWarning'].append("Annovar : transcript used is %s (instead of %s)" % (cdesc[0],annotation[variantAnnotationID]['transcript']))
					annotation[variantAnnotationID]['region'] = line['Func.refGeneWithVer']
					annotation[variantAnnotationID]['consequence'] = line['ExonicFunc.refGeneWithVer']
					annotation[variantAnnotationID]['annovar_transcriptDescription'] = cdesc[1]
					break

		# IF TRANSCRIPT NOT FOUND
		if not transcript_found:
			annotation[variantAnnotationID]['annoWarning'].append("Annovar : transcript not found")

		# COMPARE CPOS GENERATED BY HGVS MODULE VS ANNOVAR
		if transcript_found and (annotation[variantAnnotationID]['annovar_transcriptDescription'] != annotation[variantAnnotationID]['transcriptDescription']):
			annotation[variantAnnotationID]['annoWarning'].append("Annovar c. is %s" % annotation[variantAnnotationID]['annovar_transcriptDescription'])
			annotation[variantAnnotationID]['commentaire'].append("Warning : Annovar c. is %s" % annotation[variantAnnotationID]['annovar_transcriptDescription'])

		# CHANGE CONSEQUENCE DESCRIPTION :  nonsynonymous -> missense, startgain -> initiating-methionine, stopgain -> nonsense
		if annotation[variantAnnotationID]['consequence'] != None:
			if annotation[variantAnnotationID]['consequence'] == 'nonsynonymous SNV':
				annotation[variantAnnotationID]['consequence'] = 'missense'
			elif annotation[variantAnnotationID]['consequence'] == 'synonymous SNV':
				annotation[variantAnnotationID]['consequence'] = 'synonymous'
			elif 'nonframeshift' in annotation[variantAnnotationID]['consequence']:
				annotation[variantAnnotationID]['consequence'] = 'nonframeshift'
			elif 'frameshift' in annotation[variantAnnotationID]['consequence']:
				annotation[variantAnnotationID]['consequence'] = 'frameshift'

		# CORRECT REFERENCE ERROR
		if annotation[variantAnnotationID]['transcriptDescription'] != None:
			if '=' in annotation[variantAnnotationID]['transcriptDescription'] and annotation[variantAnnotationID]['consequence'] != 'synonymous':
				annotation[variantAnnotationID]['consequence'] = 'synonymous'
				annotation[variantAnnotationID]['commentaire'].append("Reference Error. Consequence switched to synonymous")

		# PARSE DATABASE ANNOTATION
		if line['cosmic92'] != '.':
			annotation[variantAnnotationID]['cosmic'] = line['cosmic92'].split('ID=')[-1].split(';OCCURENCE=')[0]
		if line['avsnp150'] != '.':
			annotation[variantAnnotationID]['dbsnp'] = line['avsnp150']
		if line['InterVar_automated'] != '.':
			annotation[variantAnnotationID]['intervar'] = line['InterVar_automated'].lower()
		if line['CLNSIG'] != '.':
			annotation[variantAnnotationID]['clinsig'] = line['CLNSIG'].lower()
		if line['nci60'] != '.':
			annotation[variantAnnotationID]['nci60'] = representsFloat(line['nci60'])
		if line['esp6500siv2_all'] != '.':
			annotation[variantAnnotationID]['esp'] = representsFloat(line['esp6500siv2_all'])
		if line['1000g2015aug_all'] != '.':
			annotation[variantAnnotationID]['1000G_all'] = representsFloat(line['1000g2015aug_all'])
		if line['1000g2015aug_eur'] != '.':
			annotation[variantAnnotationID]['1000G_eur'] = representsFloat(line['1000g2015aug_eur'])
		if line['AF'] != '.':
			annotation[variantAnnotationID]['gnomad'] = representsFloat(line['AF'])
		if line['ExAC_ALL'] != '.':
			annotation[variantAnnotationID]['exac'] = representsFloat(line['ExAC_ALL'])
		if line['SIFT_pred'] != '.':
			annotation[variantAnnotationID]['sift'] = line['SIFT_pred']  	#  D: Deleterious (sift<=0.05); T: tolerated (sift>0.05)
			if annotation[variantAnnotationID]['sift'] == 'D':
				annotation[variantAnnotationID]['sift'] = 'Deleterious'
			elif annotation[variantAnnotationID]['sift'] == 'T':
				annotation[variantAnnotationID]['sift'] = 'Tolerated'
		if line['Polyphen2_HVAR_pred'] != '.':
			annotation[variantAnnotationID]['polyphen_hvar'] = line['Polyphen2_HVAR_pred'] # "D" ("probably damaging"), "P" ("possibly damaging") and "B" ("benign").
			if annotation[variantAnnotationID]['polyphen_hvar'] == 'D':
				annotation[variantAnnotationID]['polyphen_hvar'] = 'Probably damaging'
			elif annotation[variantAnnotationID]['polyphen_hvar'] == 'P':
				annotation[variantAnnotationID]['polyphen_hvar'] = 'Possibly damaging'
			elif annotation[variantAnnotationID]['polyphen_hvar'] == 'B':
				annotation[variantAnnotationID]['polyphen_hvar'] = 'Benign'
		if line['PROVEAN_pred'] != '.':
			annotation[variantAnnotationID]['provean'] = line['PROVEAN_pred'] # "D" ("Deleterious"), "N" ("Neutral")
			if annotation[variantAnnotationID]['provean'] == 'D':
				annotation[variantAnnotationID]['provean'] = 'Deleterious'
			elif annotation[variantAnnotationID]['provean'] == 'N':
				annotation[variantAnnotationID]['provean'] = 'Neutral'

#         ___  __  
#   \  / |__  |__) 
#    \/  |___ |    

print "\n- [%s] parsing VEP results and comparing with Annovar..." % (time.strftime("%H:%M:%S"))

# PARSE VEP OUTPUT
header = []
vep_output = []

vep_transcripts = {}
wanted_transcripts = []
db_cur.execute("SELECT transcriptID from Transcript")
db_transcripts = db_cur.fetchall()
for db_transcript in db_transcripts:
	wanted_transcripts.append(db_transcript['transcriptID'])
	vep_transcripts[db_transcript['transcriptID'].split('.')[0]] = {'version':1,'searched_time':0,'stop_searching':False}

for vepline in vep_reader:
	if vepline[0].startswith('##'): # SKIP METADATA
		continue
	elif vepline[0].startswith('#'): # HEADER, USED FOR INDEX
		vepline[0].replace('#','')
		for i in range(len(vepline)):
			header.append(vepline[i])
		continue
	else:
		transcript = vepline[header.index('Feature')].split('.')[0]
		version = int(vepline[header.index('Feature')].split('.')[-1])
		if transcript in vep_transcripts.keys():
			if vep_transcripts[transcript]['stop_searching'] is True:
				continue
			else:
				if vep_transcripts[transcript]['searched_time'] >= 10 :
					vep_transcripts[transcript]['stop_searching'] = True
				elif vepline[header.index('Feature')] in wanted_transcripts: # Si dispo, prendre transcript voulu (meme si inferieur) sinon prendre "meilleur" (= le plus haut)
					vep_transcripts[transcript]['version'] = version
					vep_transcripts[transcript]['stop_searching'] = True
				elif version > vep_transcripts[transcript]['version']:
					vep_transcripts[transcript]['version'] = version
					vep_transcripts[transcript]['searched_time'] += 1

# CLEAN VEP OUTPUT (keep one line per variantAnnotation with most appropriate transcript)
vep_transcripts_to_use = ['%s.%s' % (transcript,vep_transcripts[transcript]['version']) for transcript in vep_transcripts.keys()]
vep_file.seek(0)
for vepline in vep_reader:
	if vepline[0].startswith('##'):
		continue
	elif vepline[0].startswith('#'):
		continue
	else:
		line = {}
		if '_dupl' in vepline[header.index('Feature')]: # see https://www.biostars.org/p/312593/ , bug duplicate
			continue
		if not vepline[header.index('Feature')].startswith('NM_'): # Exclude NR, XM ...
			continue
		if vepline[header.index('Feature')] in vep_transcripts_to_use:
			for i in range(len(vepline)):
				line[header[i]] = vepline[i]
			vep_output.append(line)

# PARCOURS DES VARIANTS VEP
z=0
for vepline in vep_output:
	z+=1
	print "\r\t - %s/%s" % (z,len(vep_output)),
	vep_transcript = vepline['Feature']
	vep_transcript_without_version = vepline['Feature'].split('.')[0]

	# DETERMINE VARIANT ID
	location = vepline['Location']
	chrom = 'chr%s' % location.split(':')[0]
	ref = vepline['GIVEN_REF']
	alt = vepline['Allele']
	if '-' in location :
		start = location.split(':')[-1].split('-')[0]
		stop = location.split(':')[-1].split('-')[-1]
	else:
		start = location.split(':')[-1]
		stop = start
	variantID = '%s:%s-%s:%s>%s' % (chrom,start,stop,ref,alt)

	db_cur.execute("SELECT * FROM VariantAnnotation INNER JOIN Variant ON Variant.variantID=VariantAnnotation.variant WHERE variant='%s'" % variantID)
	db_variantannotations = db_cur.fetchall() # il peut y avoir plusieurs entrees si plusieurs transcrits pour un meme gene
	for db_variantannotation in db_variantannotations:
		variantAnnotationID = db_variantannotation['variantAnnotationID']
		if variantAnnotationID not in annotation.keys():
			print "\n\t-VEP Warning : VEP variant %s does not match any ANNOVAR variant" % vepline['Uploaded_variation']
			continue

		if annotation[variantAnnotationID]['transcript_without_version'] == vep_transcript_without_version :
			if annotation[variantAnnotationID]['transcript'] == vep_transcript :
				annotation[variantAnnotationID]['annoWarning'].append("VEP : transcript used is %s (instead of %s)" % (vep_transcript,annotation[variantAnnotationID]['transcript']))
			if vepline['EXON'] != '-':
				annotation[variantAnnotationID]['vep_exon'] = representsInt(vepline['EXON'].split('/')[0])
			if vepline['INTRON'] != '-':
				annotation[variantAnnotationID]['vep_intron'] = representsInt(vepline['INTRON'].split('/')[0])
			if vepline['HGVSc'] != '-':
				annotation[variantAnnotationID]['vep_transcriptDescription'] = vepline['HGVSc'].split(':')[-1]
				if annotation[variantAnnotationID]['vep_transcriptDescription'] != annotation[variantAnnotationID]['transcriptDescription']:
					annotation[variantAnnotationID]['annoWarning'].append("VEP c. is %s" % annotation[variantAnnotationID]['vep_transcriptDescription'])
			if vepline['Consequence'] != '-':
				annotation[variantAnnotationID]['vep_consequence'] = vepline['Consequence']
			if vepline['IMPACT'] != '-':
				annotation[variantAnnotationID]['vep_impact'] = vepline['IMPACT']
			if vepline['PUBMED'] != '-':
				annotation[variantAnnotationID]['pubmed'] = vepline['PUBMED']
			vep_dbSNP = []
			vep_COSMIC = []
			if vepline['Existing_variation'] != '-':
				existing_variation = vepline['Existing_variation'].split(',')
				for item in existing_variation:
					if item.startswith('rs'):
						vep_dbSNP.append(item)
					elif item.startswith('COSM'):
						vep_COSMIC.append(item)
				if vep_dbSNP:
					annotation[variantAnnotationID]['vep_dbsnp'] = ','.join(vep_dbSNP)
				if vep_COSMIC:
					annotation[variantAnnotationID]['vep_cosmic'] = ','.join(vep_COSMIC)
			if vepline['CLIN_SIG'] != '-':
				annotation[variantAnnotationID]['vep_clinsig'] = vepline['CLIN_SIG']
			if vepline['EA_AF'] != '-' and vepline['AA_AF'] != '-':
				try:
					annotation[variantAnnotationID]['vep_esp'] = str((float(vepline['EA_AF']) + float(vepline['AA_AF']))/2)
				except:
					annotation[variantAnnotationID]['vep_esp'] = None
			if vepline['AF'] != '-':
				annotation[variantAnnotationID]['vep_1000G_all'] = representsFloat(vepline['AF'])
			if vepline['EUR_AF'] != '-':
				annotation[variantAnnotationID]['vep_1000G_eur'] = representsFloat(vepline['EUR_AF'])
			if vepline['SIFT'] != '-':
				annotation[variantAnnotationID]['vep_sift'] = vepline['SIFT']
			if vepline['PolyPhen'] != '-':
				annotation[variantAnnotationID]['vep_polyphen_hvar'] = vepline['PolyPhen']

			#########################################
			## DETERMINING ANNOVAR/VEP DIFFERENCES ##		## N'est pas present dans VEP : INTERVAR, NCI60, GNOMAD, EXAC, PROVEAN
			#########################################

			diff = {
				'exon':[],
				'clinsig':[],
				'cosmic':[],
				'dbsnp':[],
				'1000G_ALL':[],
				'1000G_EUR':[],
				'esp':[],
				'sift':[],
				'polyphen_hvar':[]
			}

			## EXON DIFF (une valeur possible)
			if annotation[variantAnnotationID]['exon'] != None and annotation[variantAnnotationID]['vep_exon'] != None:
				if annotation[variantAnnotationID]['exon'] != annotation[variantAnnotationID]['vep_exon']:
					diff['exon'].append('+%s'%annotation[variantAnnotationID]['vep_exon'])
					diff['exon'].append('-%s'%annotation[variantAnnotationID]['exon'])
			elif annotation[variantAnnotationID]['exon'] == None and annotation[variantAnnotationID]['vep_exon'] != None:
				diff['exon'].append('+%s'%annotation[variantAnnotationID]['vep_exon'])
			elif annotation[variantAnnotationID]['exon'] != None and annotation[variantAnnotationID]['vep_exon'] == None:
				diff['exon'].append('-%s'%annotation[variantAnnotationID]['exon'])

			## CLINVAR DIFF (plusieurs valeurs possibles)
			annovar_clinsig = []
			vep_clinsig = []
			if annotation[variantAnnotationID]['clinsig'] != None:
				annovar_clinsig = annotation[variantAnnotationID]['clinsig'].replace('/',',').split(',')
			if annotation[variantAnnotationID]['vep_clinsig'] != None:
				vep_clinsig = annotation[variantAnnotationID]['vep_clinsig'].split(',')
			if 'not_provided' in vep_clinsig:
				vep_clinsig.remove('not_provided')
			for cs in vep_clinsig:
				if cs not in annovar_clinsig: # and cs != '-'
					diff['clinsig'].append('+%s'%cs)
			for cs in annovar_clinsig:
				if cs not in vep_clinsig: #  and cs != '.'
					diff['clinsig'].append('-%s'%cs)

			if annotation[variantAnnotationID]['clinsig'] == None and annotation[variantAnnotationID]['vep_clinsig'] != None: ## Si Annovar == None, utilisation VEP par defaut
				if 'benign' in annotation[variantAnnotationID]['vep_clinsig'] and 'pathogenic' in annotation[variantAnnotationID]['vep_clinsig']:
					annotation[variantAnnotationID]['clinsig'] = 'Conflicting_interpretations_of_pathogenicity'
				else:
					if 'pathogenic' in vep_clinsig and 'likely_pathogenic' in vep_clinsig:
						vep_clinsig.remove('pathogenic')
						vep_clinsig.remove('likely_pathogenic')
						vep_clinsig.insert(0,'pathogenic/likely_pathogenic')
					if 'likely_pathogenic' in vep_clinsig:
						vep_clinsig.remove('likely_pathogenic')
						vep_clinsig.insert(0,'likely_pathogenic')
					if 'pathogenic' in vep_clinsig:
						vep_clinsig.remove('pathogenic')
						vep_clinsig.insert(0,'pathogenic')
					annotation[variantAnnotationID]['clinsig'] = ','.join(vep_clinsig)

			## COSMIC DIFF (plusieurs valeurs possibles)
			annovar_COSMIC = []
			if annotation[variantAnnotationID]['cosmic'] != None:
				annovar_COSMIC = annotation[variantAnnotationID]['cosmic'].split(',')
			for cosm in vep_COSMIC:
				if cosm not in annovar_COSMIC: #  and cosm != '.'
					diff['cosmic'].append('+%s'%cosm)
			for cosm in annovar_COSMIC:
				if cosm not in vep_COSMIC: #  and cosm != '.'
					diff['cosmic'].append('-%s'%cosm)

			## DBSNP DIFF (une valeur possible)
			if annotation[variantAnnotationID]['dbsnp'] != None and annotation[variantAnnotationID]['vep_dbsnp'] != None:
				if annotation[variantAnnotationID]['dbsnp'] != annotation[variantAnnotationID]['vep_dbsnp']:
					diff['dbsnp'].append('+%s'%annotation[variantAnnotationID]['vep_dbsnp'])
					diff['dbsnp'].append('-%s'%annotation[variantAnnotationID]['dbsnp'])
			elif annotation[variantAnnotationID]['dbsnp'] == None and annotation[variantAnnotationID]['vep_dbsnp'] != None:
				diff['dbsnp'].append('+%s'%annotation[variantAnnotationID]['vep_dbsnp'])
			elif annotation[variantAnnotationID]['dbsnp'] != None and annotation[variantAnnotationID]['vep_dbsnp'] == None:
				diff['dbsnp'].append('-%s'%annotation[variantAnnotationID]['dbsnp'])

			## 1000G_ALL DIFF (une valeur possible) # 2% de tolerance
			if annotation[variantAnnotationID]['1000G_all'] != None and annotation[variantAnnotationID]['vep_1000G_all'] != None:
				if not ((round(float(annotation[variantAnnotationID]['1000G_all']),2)-0.02) <= round(float(annotation[variantAnnotationID]['vep_1000G_all']),2) <= (round(float(annotation[variantAnnotationID]['1000G_all']),2)+0.02)):
					diff['1000G_ALL'].append('+%s'%annotation[variantAnnotationID]['vep_1000G_all'])
					diff['1000G_ALL'].append('-%s'%annotation[variantAnnotationID]['1000G_all'])
			elif annotation[variantAnnotationID]['1000G_all'] == None and annotation[variantAnnotationID]['vep_1000G_all'] != None:
				diff['1000G_ALL'].append('+%s'%annotation[variantAnnotationID]['vep_1000G_all'])
			elif annotation[variantAnnotationID]['1000G_all'] != None and annotation[variantAnnotationID]['vep_1000G_all'] == None:
				diff['1000G_ALL'].append('-%s'%annotation[variantAnnotationID]['1000G_all'])

			## 1000G_EUR DIFF (une valeur possible) # 2% de tolerance
			if annotation[variantAnnotationID]['1000G_eur'] != None and annotation[variantAnnotationID]['vep_1000G_eur'] != None:
				if not ((round(float(annotation[variantAnnotationID]['1000G_eur']),2)-0.02) <= round(float(annotation[variantAnnotationID]['vep_1000G_eur']),2) <= (round(float(annotation[variantAnnotationID]['1000G_eur']),2)+0.02)):
					diff['1000G_EUR'].append('+%s'%annotation[variantAnnotationID]['vep_1000G_eur'])
					diff['1000G_EUR'].append('-%s'%annotation[variantAnnotationID]['1000G_eur'])
			elif annotation[variantAnnotationID]['1000G_eur'] == None and annotation[variantAnnotationID]['vep_1000G_eur'] != None:
				diff['1000G_EUR'].append('+%s'%annotation[variantAnnotationID]['vep_1000G_eur'])
			elif annotation[variantAnnotationID]['1000G_eur'] != None and annotation[variantAnnotationID]['vep_1000G_eur'] == None:
				diff['1000G_EUR'].append('-%s'%annotation[variantAnnotationID]['1000G_eur'])

			## ESP DIFF (une valeur possible) # 2% de tolerance
			if annotation[variantAnnotationID]['esp'] != None and annotation[variantAnnotationID]['vep_esp'] != None:	
				if not ((round(float(annotation[variantAnnotationID]['esp']),2)-0.02) <= round(float(annotation[variantAnnotationID]['vep_esp']),2) <= (round(float(annotation[variantAnnotationID]['esp']),2)+0.02)):
					diff['esp'].append('+%s'%annotation[variantAnnotationID]['vep_esp'])
					diff['esp'].append('-%s'%annotation[variantAnnotationID]['esp'])
			elif annotation[variantAnnotationID]['esp'] == None and annotation[variantAnnotationID]['vep_esp'] != None:
				diff['esp'].append('+%s'%annotation[variantAnnotationID]['vep_esp'])
			elif annotation[variantAnnotationID]['esp'] != None and annotation[variantAnnotationID]['vep_esp'] == None:
				diff['esp'].append('-%s'%annotation[variantAnnotationID]['esp'])

			## SIFT DIFF (une valeur possible)
			if annotation[variantAnnotationID]['sift'] != None and annotation[variantAnnotationID]['vep_sift'] != None:	
				if str(annotation[variantAnnotationID]['sift']).lower() != str(annotation[variantAnnotationID]['vep_sift']).replace('_',' ').replace('-','.'):
					diff['sift'].append('+%s'%annotation[variantAnnotationID]['vep_sift'])
					diff['sift'].append('-%s'%annotation[variantAnnotationID]['sift'])
			elif annotation[variantAnnotationID]['sift'] == None and annotation[variantAnnotationID]['vep_sift'] != None:
				diff['sift'].append('+%s'%annotation[variantAnnotationID]['vep_sift'])
			elif annotation[variantAnnotationID]['sift'] != None and annotation[variantAnnotationID]['vep_sift'] == None:
				diff['sift'].append('-%s'%annotation[variantAnnotationID]['sift'])

			## POLYPHEN DIFF (une valeur possible)
			if annotation[variantAnnotationID]['polyphen_hvar'] != None and annotation[variantAnnotationID]['vep_polyphen_hvar'] != None:	
				if str(annotation[variantAnnotationID]['polyphen_hvar']).lower() != str(annotation[variantAnnotationID]['vep_polyphen_hvar']).replace('_',' ').replace('-','.'):
					diff['polyphen_hvar'].append('+%s'%annotation[variantAnnotationID]['vep_polyphen_hvar'])
					diff['polyphen_hvar'].append('-%s'%annotation[variantAnnotationID]['polyphen_hvar'])
			elif annotation[variantAnnotationID]['polyphen_hvar'] == None and annotation[variantAnnotationID]['vep_polyphen_hvar'] != None:
				diff['polyphen_hvar'].append('+%s'%annotation[variantAnnotationID]['vep_polyphen_hvar'])
			elif annotation[variantAnnotationID]['polyphen_hvar'] != None and annotation[variantAnnotationID]['vep_polyphen_hvar'] == None:
				diff['polyphen_hvar'].append('-%s'%annotation[variantAnnotationID]['polyphen_hvar'])

			#### GENERATE DIFF STRING ####
			vep_diff = []
			for item in diff.keys():
				if diff[item]:
					d = ','.join(diff[item])
					vep_diff.append('%s:%s'%(item,d))
			if vep_diff:
				 annotation[variantAnnotationID]['vep_diff'] = ';'.join(vep_diff)#annotation[variantAnnotationID]['vep_diff'][1:]

########################################
## COMPARING ANNOVAR/VEP CONSEQUENCES ##
########################################

print "\n- [%s] comparing Annovar/VEP consequences ..." % (time.strftime("%H:%M:%S"))
z=0
for variantAnnotationID in annotation.keys():
	z+=1
	print "\r\t - %s/%s" % (z,len(annotation)),
	if annotation[variantAnnotationID]['consequence'] == None:
		if annotation[variantAnnotationID]['vep_consequence'] != None:
			vep_consequences = annotation[variantAnnotationID]['vep_consequence'].split(',') # sometimes several consequences
			region_consequence_list = []
			for consequence in vep_consequences:
				if consequence in vep2annovar:
					region_consequence_list.append((vep2annovar[consequence]['region'],vep2annovar[consequence]['type']))
			final_consequence = compare_consequence(region_consequence_list)
			annotation[variantAnnotationID]['region'] = final_consequence[0]
			annotation[variantAnnotationID]['consequence'] = final_consequence[1]
			annotation[variantAnnotationID]['dbsnp'] = annotation[variantAnnotationID]['vep_dbsnp']
			annotation[variantAnnotationID]['cosmic'] = annotation[variantAnnotationID]['vep_cosmic']
			annotation[variantAnnotationID]['dbsnp'] = annotation[variantAnnotationID]['vep_dbsnp']
			annotation[variantAnnotationID]['clinsig'] = annotation[variantAnnotationID]['vep_clinsig']
			annotation[variantAnnotationID]['esp'] = annotation[variantAnnotationID]['vep_esp']
			annotation[variantAnnotationID]['1000G_all'] = annotation[variantAnnotationID]['vep_1000G_all']
			annotation[variantAnnotationID]['sift'] = annotation[variantAnnotationID]['vep_sift']
			annotation[variantAnnotationID]['polyphen_hvar'] = annotation[variantAnnotationID]['vep_polyphen_hvar']
			annotation[variantAnnotationID]['annoWarning'].append("VEP annotation was used")
	elif annotation[variantAnnotationID]['vep_consequence'] != None:
		vep_consequences = annotation[variantAnnotationID]['vep_consequence'].split(',') # sometimes several consequences
		region_consequence_list = [(annotation[variantAnnotationID]['region'],annotation[variantAnnotationID]['consequence'])]
		for consequence in vep_consequences:
			if consequence in vep2annovar:
				region_consequence_list.append((vep2annovar[consequence]['region'],vep2annovar[consequence]['type']))
		final_consequence = compare_consequence(region_consequence_list)
		if final_consequence != (annotation[variantAnnotationID]['region'],annotation[variantAnnotationID]['consequence']):
			annotation[variantAnnotationID]['annoWarning'].append("Annovar/VEP discordance for consequence")
	# TODO : modifier vep2annovar et enlever les '.', en attendant:
	if annotation[variantAnnotationID]['consequence'] == '.':
		annotation[variantAnnotationID]['consequence'] = None


#############################
# ADDING CUSTOM ANNOTATIONS #
#############################

# KNOWN FALSE-POSITIVES, DRUG SENSITIVITY, LAM HOTSPOT HIGHLIGHTS, TP53 UMD
# print "\n- [%s] adding custom annotations ..." % (time.strftime("%H:%M:%S"))

# custom_fp = []
# with open(custom_false_positives,'r') as fp:
	# for row in fp :
		# fp_row = row.replace('\n','').split(';')
		# custom_fp.append(fp_row)

# sensi = []
# sensitivity_xls = xlrd.open_workbook(custom_drug_sensitivity)
# sensitivity_sheet = sensitivity_xls.sheet_by_index(0)
# for row in range(sensitivity_sheet.nrows-1):
	# sensi.append((str(sensitivity_sheet.cell(row+1,5).value),str(sensitivity_sheet.cell(row+1,7).value),str(sensitivity_sheet.cell(row+1,8).value),str(sensitivity_sheet.cell(row+1,9).value))) # gene,c,p,sensi

# custom_hl = []
# with open(custom_highlight,'r') as hl:
	# for row in hl:
		# hl_row = row.replace('\n','')
		# hl_row = hl_row.split(';')
		# custom_hl.append((hl_row[0].split('.')[0],hl_row[1])) # NM without version

# custom_tu = []
# with open(custom_tp53_umd,'r') as tp53_umd:
	# tp53_umd.next()
	# for row in tp53_umd:
		# tp53_row = row.replace('\n','')
		# tp53_row = tp53_row.split('\t')
		# custom_tu.append((tp53_row[0],tp53_row[1],tp53_row[2],tp53_row[3]))

# custom_lc = []
# with open(custom_lymphome_controls,'r') as lc:
	# for row in lc :
		# lc_row = row.replace('\n','').split('\t')
		# nb = len(lc_row[2].split(','))
		# string = "Found in normal dna (x%s)" % nb # "Found in %s controls (%s)" % (nb,lc_row[2])
		# lc_row[2] = string
		# custom_lc.append(lc_row)

# z=0
# for variantAnnotationID in annotation.keys():
	# z+=1
	# print "\r\t - %s/%s" % (z,len(annotation)),
	# # KNOWN FALSE-POSITIVES
	# for fp in custom_fp:
		# if (fp[0] == annotation[variantAnnotationID]['transcript_without_version']) and ((fp[1] == annotation[variantAnnotationID]['transcriptDescription']) or (fp[1] == annotation[variantAnnotationID]['annovar_transcriptDescription'])):
			# # annotation[variantAnnotationID]['commentaire'].append(fp[2])
			# annotation[variantAnnotationID]['commentaire'].insert(0,fp[2])
			# break
	# # DRUG SENSITIVITY
	# for s in sensi:
		# if s[2] not in ['','p.(=)','p.?','p.=']:
			# if ((annotation[variantAnnotationID]['gene'],annotation[variantAnnotationID]['proteinDescription']) == (s[0],s[2])) or ((annotation[variantAnnotationID]['gene'],annotation[variantAnnotationID]['annovar_proteinDescription']) == (s[0],s[2])) :
				# annotation[variantAnnotationID]['actionability'] = s[3]
				# break
		# if s[1] != '':
			# if ((annotation[variantAnnotationID]['gene'],annotation[variantAnnotationID]['transcriptDescription']) == (s[0],s[1])) or ((annotation[variantAnnotationID]['gene'],annotation[variantAnnotationID]['annovar_transcriptDescription']) == (s[0],s[1])):
				# annotation[variantAnnotationID]['actionability'] = s[3]
				# break
	# ## cas particulier : splicing MET
	# if (annotation[variantAnnotationID]['gene'] == 'MET') and (annotation[variantAnnotationID]['region'] == 'intronic') and ((116411873 <= annotation[variantAnnotationID]['genomicStart'] <= 116411902) or (116412044 <= annotation[variantAnnotationID]['genomicStart'] <= 116412087)):
		# annotation[variantAnnotationID]['exon'] = 14
		# annotation[variantAnnotationID]['commentaire'].append('MET intron 13/14')
		# if annotation[variantAnnotationID]['actionability'] == None :
			# annotation[variantAnnotationID]['actionability'] = 'intron 13/14'
	# # LAM HOTSPOT HIGHLIGHTS
	# for hl in custom_hl:
		# if ((annotation[variantAnnotationID]['transcript_without_version'],annotation[variantAnnotationID]['transcriptDescription']) == hl) or ((annotation[variantAnnotationID]['transcript_without_version'],annotation[variantAnnotationID]['annovar_transcriptDescription']) == hl): # or ((annotation[variantAnnotationID]['gene'],annotation[variantAnnotationID]['transcriptDescription'],annotation[variantAnnotationID]['proteinDescription']) == hl) # necessite hgvs avant
			# annotation[variantAnnotationID]['highlight'] = 1
			# break
	# #TP53 UMD
	# for tu in custom_tu:
		# if tu[1] not in ['','p.(=)','p.?','p.=']:
			# if ((annotation[variantAnnotationID]['gene'],annotation[variantAnnotationID]['proteinDescription']) == ('TP53',tu[1])) or ((annotation[variantAnnotationID]['gene'],annotation[variantAnnotationID]['annovar_proteinDescription']) == ('TP53',tu[1])):
				# annotation[variantAnnotationID]['pathoUMD'] = tu[2]
				# annotation[variantAnnotationID]['commentUMD'] = tu[3]
				# break
		# if tu[0] not in ['','c.(=)','c.?','c.=']:
			# if ((annotation[variantAnnotationID]['gene'],annotation[variantAnnotationID]['transcriptDescription']) == ('TP53',tu[0])) or ((annotation[variantAnnotationID]['gene'],annotation[variantAnnotationID]['annovar_transcriptDescription']) == ('TP53',tu[0])):
				# annotation[variantAnnotationID]['pathoUMD'] = tu[2]
				# annotation[variantAnnotationID]['commentUMD'] = tu[3]
				# break
	# #LYMPHOME CONTROLS
	# for lc in custom_lc:
		# if (lc[0] == annotation[variantAnnotationID]['transcript_without_version']) and ((lc[1] == annotation[variantAnnotationID]['transcriptDescription']) or (lc[1] == annotation[variantAnnotationID]['annovar_transcriptDescription'])):
			# annotation[variantAnnotationID]['commentaire'].append(lc[2])
			# break

#################
# WRITE RESULTS #
#################

print "\n- [%s] writting results to VariantBase ..." % (time.strftime("%H:%M:%S"))
lastUpdate = '%s' % (time.strftime('%d/%m/%Y'))
for variantAnnotationID in annotation.keys():
	if annotation[variantAnnotationID]['commentaire']:
		commentaire = '. '.join(annotation[variantAnnotationID]['commentaire'])
	else:
		commentaire = None
	if annotation[variantAnnotationID]['annoWarning']:
		annoWarning = '. '.join(annotation[variantAnnotationID]['annoWarning'])
	else:
		annoWarning = None

	try:
		db_cur.execute("""UPDATE VariantAnnotation SET exon=?,intron=?,region=?,consequence=?,actionability=?,intervar=?,clinvar=?,cosmic=?,dbsnp=?,gnomad=?,milleGall=?,milleGeur=?,nci60=?,esp=?,exac=?,sift=?,polyphen2=?,provean=?,pubmed=?,vep_consequence=?,vep_impact=?,vep_diff=?,commentaire=?,highlight=?,pathoUMD=?,commentUMD=?,annovarTranscriptDescription=?,annovarProteinDescription=?,annoWarning=?,lastUpdate=? WHERE variantAnnotationID=?""", (annotation[variantAnnotationID]['exon'],annotation[variantAnnotationID]['vep_intron'],annotation[variantAnnotationID]['region'],annotation[variantAnnotationID]['consequence'],annotation[variantAnnotationID]['actionability'],annotation[variantAnnotationID]['intervar'],annotation[variantAnnotationID]['clinsig'],annotation[variantAnnotationID]['cosmic'],annotation[variantAnnotationID]['dbsnp'],annotation[variantAnnotationID]['gnomad'],annotation[variantAnnotationID]['1000G_all'],annotation[variantAnnotationID]['1000G_eur'],annotation[variantAnnotationID]['nci60'],annotation[variantAnnotationID]['esp'],annotation[variantAnnotationID]['exac'],annotation[variantAnnotationID]['sift'],annotation[variantAnnotationID]['polyphen_hvar'],annotation[variantAnnotationID]['provean'],annotation[variantAnnotationID]['pubmed'],annotation[variantAnnotationID]['vep_consequence'],annotation[variantAnnotationID]['vep_impact'],annotation[variantAnnotationID]['vep_diff'],commentaire,annotation[variantAnnotationID]['highlight'],annotation[variantAnnotationID]['pathoUMD'],annotation[variantAnnotationID]['commentUMD'],annotation[variantAnnotationID]['annovar_transcriptDescription'],annotation[variantAnnotationID]['annovar_proteinDescription'],annoWarning,lastUpdate,variantAnnotationID))
	except Exception as e:
		print "Error : %s" % e

db_con.commit()
db_con.close()
print "- [%s] ... done." % (time.strftime("%H:%M:%S"))
