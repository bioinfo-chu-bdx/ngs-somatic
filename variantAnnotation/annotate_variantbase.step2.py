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
	print "- 0 new variants to process"
	exit()

annovar_file = open(options.annovar_results,'r')
annovar_reader = csv.DictReader(annovar_file,delimiter='\t')

vep_file = open(options.vep_results,'r')
vep_reader = csv.reader(vep_file,delimiter='\t')

custom_false_positives = global_param['fp_file']
custom_drug_sensitivity = global_param['run_type']['SBT']['sbt_sensitivity']
custom_highlight = global_param['run_type']['LAM']['lam_hotspot_highlight']
custom_tp53_umd = global_param['run_type']['TP53']['TP53_UMD_variants_EU']
custom_lymphome_controls = global_param['run_type']['Lymphome_T']['lymphomeControls']

annotation = {}

gene_data = {}
db_cur.execute("SELECT * FROM Gene")
db_genes = db_cur.fetchall()
for g in db_genes:
	gene_data[g['geneID']] = {'transcript': g['transcript'],'transcriptVersion': g['transcriptVersion']}
	
############################
# PARSE ANNOTATION RESULTS #
############################
#                   __             __  
#    /\  |\ | |\ | /  \ \  /  /\  |__) 
#   /~~\ | \| | \| \__/  \/  /~~\ |  \ 

print "- [%s] parsing annovar results ..." % (time.strftime("%H:%M:%S"))
for line in annovar_reader:
	variant = '%s:%s-%s:%s>%s' % (line['Chr'],line['Start'],line['End'],line['Ref'],line['Alt'])
	if line['Ref'] == '-':
		variant = '%s:%s-%s:%s>%s' % (line['Chr'],line['Start'],str(int(line['Start'])+1),line['Ref'],line['Alt'])
		
	annotation[variant] = {}
	annotation[variant]['commentaire'] = []
	annotation[variant]['annoWarning'] = []
	annotation[variant]['gene'] = None
	annotation[variant]['genomicStart'] = None
	annotation[variant]['transcript'] = None
	annotation[variant]['transcriptVersion'] = None
	annotation[variant]['exon'] = None
	annotation[variant]['vep_exon'] = None
	annotation[variant]['vep_intron'] = None
	annotation[variant]['region'] = None
	annotation[variant]['consequence'] = None
	annotation[variant]['vep_consequence'] = None
	annotation[variant]['vep_impact'] = None
	annotation[variant]['transcriptDescription'] = None
	annotation[variant]['annovar_transcriptDescription'] = None
	annotation[variant]['vep_transcriptDescription'] = None
	annotation[variant]['proteinDescription'] = None
	annotation[variant]['annovar_proteinDescription'] = None
	annotation[variant]['actionability'] = None
	annotation[variant]['pubmed'] = None
	annotation[variant]['cosmic'] = None
	annotation[variant]['vep_cosmic'] = None
	annotation[variant]['dbsnp'] = None
	annotation[variant]['vep_dbsnp'] = None
	annotation[variant]['intervar'] = None
	annotation[variant]['clinsig'] = None
	annotation[variant]['vep_clinsig'] = None
	annotation[variant]['nci60'] = None
	annotation[variant]['esp'] = None
	annotation[variant]['vep_esp'] = None
	annotation[variant]['1000G_all'] = None
	annotation[variant]['vep_1000G_all'] = None
	annotation[variant]['1000G_eur'] = None
	annotation[variant]['vep_1000G_eur'] = None
	annotation[variant]['gnomad'] = None
	annotation[variant]['exac'] = None
	annotation[variant]['sift'] = None
	annotation[variant]['vep_sift'] = None
	annotation[variant]['polyphen_hvar'] = None
	annotation[variant]['vep_polyphen_hvar'] = None
	annotation[variant]['provean'] = None
	annotation[variant]['highlight'] = 0
	annotation[variant]['pathoUMD'] = None
	annotation[variant]['commentUMD'] = None
	annotation[variant]['vep_diff'] = None
	
	db_cur.execute("SELECT * FROM Variant WHERE variantID='%s'" % variant)
	db_variant = db_cur.fetchone()
	gene = db_variant['gene']
	annotation[variant]['gene'] = gene
	annotation[variant]['transcript'] = gene_data[gene]['transcript']
	annotation[variant]['transcriptVersion'] = gene_data[gene]['transcriptVersion']
	transcrit_with_version = '%s.%s' % (gene_data[gene]['transcript'],gene_data[gene]['transcriptVersion'])
	annotation[variant]['transcriptDescription'] = db_variant['transcriptDescription']
	annotation[variant]['proteinDescription'] = db_variant['proteinDescription']
	annotation[variant]['genomicStart'] = db_variant['genomicStart']
	
	nm_found = False
	if not line['AAChange.refGeneWithVer'] == '.':
		aachange = line['AAChange.refGeneWithVer'].split(',') # format GENE:NM:exon:c:p,
		for item in aachange:
			aadesc = item.split(':')
			if annotation[variant]['transcript'] in aadesc[1]:
				if aadesc[1] != transcrit_with_version:
					annotation[variant]['annoWarning'].append("Annovar transcript version : '%s'" % aadesc[1])
					#print "\t-Annovar transcript version is '%s'" % aadesc[1]
				nm_found = True
				annotation[variant]['region'] = line['Func.refGeneWithVer']
				annotation[variant]['consequence'] = line['ExonicFunc.refGeneWithVer']
				annotation[variant]['exon'] = representsInt(aadesc[2].replace('exon',''))
				annotation[variant]['annovar_transcriptDescription'] = aadesc[3]
				annotation[variant]['annovar_proteinDescription'] = aadesc[4]
				break
	# non-exonic avec c. (cas particulier : aachange est vide et l'equivalent se trouve dans dans genedetail)
	if line['GeneDetail.refGeneWithVer'] != '.':
		if 'splicing' in line['Func.refGeneWithVer']:
			cchange = line['GeneDetail.refGeneWithVer'].split(';') # format NM:exon:c;
			for item in cchange:
				cdesc = item.split(':')
				if annotation[variant]['transcript'] in cdesc[0]:
					nm_found = True
					annotation[variant]['region'] = line['Func.refGeneWithVer']
					annotation[variant]['consequence'] = line['ExonicFunc.refGeneWithVer']
					annotation[variant]['exon'] = representsInt(cdesc[1].replace('exon',''))
					annotation[variant]['annovar_transcriptDescription'] = cdesc[2]
					break
		else : # ncRNA, UTR3, UTR5, upstream, downstream
			cchange = line['GeneDetail.refGeneWithVer'].split(';')  # format NM:c;
			for item in cchange:
				cdesc = item.split(':') 
				if annotation[variant]['transcript'] in cdesc[0]:
					nm_found = True
					annotation[variant]['region'] = line['Func.refGeneWithVer']
					annotation[variant]['consequence'] = line['ExonicFunc.refGeneWithVer']
					annotation[variant]['annovar_transcriptDescription'] = cdesc[1]
					break
	# COMPARE c
	if nm_found and (annotation[variant]['annovar_transcriptDescription'] != annotation[variant]['transcriptDescription']):
		annotation[variant]['annoWarning'].append("Annovar : %s" % annotation[variant]['annovar_transcriptDescription'])
		#print "\t-Annovar transcriptDescription is %s" % annotation[variant]['annovar_transcriptDescription']
		
	# consequence :  nonsynonymous -> missense, startgain -> initiating-methionine, stopgain -> nonsense
	if annotation[variant]['consequence'] != None:
		if annotation[variant]['consequence'] == 'nonsynonymous SNV':
			annotation[variant]['consequence'] = 'missense'
		elif annotation[variant]['consequence'] == 'synonymous SNV':
			annotation[variant]['consequence'] = 'synonymous'
		elif 'nonframeshift' in annotation[variant]['consequence']:
			annotation[variant]['consequence'] = 'nonframeshift'
		elif 'frameshift' in annotation[variant]['consequence']:
			annotation[variant]['consequence'] = 'frameshift'
			
	## CORRECT REFERENCE ERROR
	if annotation[variant]['transcriptDescription'] != None:
		if '=' in annotation[variant]['transcriptDescription'] and annotation[variant]['consequence'] != 'synonymous':
			annotation[variant]['consequence'] = 'synonymous'
			annotation[variant]['commentaire'].append("Reference Error. Consequence switched to synonymous")
			#print "\t-Reference Error : consequence changed to synonymous"
	
	if not nm_found:
		annotation[variant]['annoWarning'].append("Annovar : no data for %s" % annotation[variant]['transcript'])
		#print "\t-Annovar Warning : Annovar could not find data for %s transcript" % annotation[variant]['transcript']

	if line['cosmic92'] != '.':
		annotation[variant]['cosmic'] = line['cosmic92'].split('ID=')[-1].split(';OCCURENCE=')[0]
	if line['avsnp150'] != '.':
		annotation[variant]['dbsnp'] = line['avsnp150']
	if line['InterVar_automated'] != '.':
		annotation[variant]['intervar'] = line['InterVar_automated'].lower()
	if line['CLNSIG'] != '.':
		annotation[variant]['clinsig'] = line['CLNSIG'].lower()
	if line['nci60'] != '.':
		annotation[variant]['nci60'] = representsFloat(line['nci60'])
	if line['esp6500siv2_all'] != '.':
		annotation[variant]['esp'] = representsFloat(line['esp6500siv2_all'])
	if line['1000g2015aug_all'] != '.':
		annotation[variant]['1000G_all'] = representsFloat(line['1000g2015aug_all'])
	if line['1000g2015aug_eur'] != '.':
		annotation[variant]['1000G_eur'] = representsFloat(line['1000g2015aug_eur'])
	if line['AF'] != '.':
		annotation[variant]['gnomad'] = representsFloat(line['AF'])
	if line['ExAC_ALL'] != '.':
		annotation[variant]['exac'] = representsFloat(line['ExAC_ALL'])

	if line['SIFT_pred'] != '.':
		annotation[variant]['sift'] = line['SIFT_pred']  	#  D: Deleterious (sift<=0.05); T: tolerated (sift>0.05)
		if annotation[variant]['sift'] == 'D':
			annotation[variant]['sift'] = 'Deleterious'
		elif annotation[variant]['sift'] == 'T':
			annotation[variant]['sift'] = 'Tolerated'
	
	if line['Polyphen2_HVAR_pred'] != '.':
		annotation[variant]['polyphen_hvar'] = line['Polyphen2_HVAR_pred'] # "D" ("probably damaging"), "P" ("possibly damaging") and "B" ("benign").
		if annotation[variant]['polyphen_hvar'] == 'D':
			annotation[variant]['polyphen_hvar'] = 'Probably damaging'
		elif annotation[variant]['polyphen_hvar'] == 'P':
			annotation[variant]['polyphen_hvar'] = 'Possibly damaging'
		elif annotation[variant]['polyphen_hvar'] == 'B':
			annotation[variant]['polyphen_hvar'] = 'Benign'
			
	if line['PROVEAN_pred'] != '.':
		annotation[variant]['provean'] = line['PROVEAN_pred'] # "D" ("Deleterious"), "N" ("Neutral")
		if annotation[variant]['provean'] == 'D':
			annotation[variant]['provean'] = 'Deleterious'
		elif annotation[variant]['provean'] == 'N':
			annotation[variant]['provean'] = 'Neutral'

#         ___  __  
#   \  / |__  |__) 
#    \/  |___ |    

print "- [%s] parsing VEP results and comparing with Annovar..." % (time.strftime("%H:%M:%S"))

vep_transcript_version = {} # in VEP results, select higher transcript version to use
index = {}
for gene in gene_data:
	vep_transcript_version[gene_data[gene]['transcript']] = 1
for vepline in vep_reader:
	if vepline[0].startswith('##'):
		continue
	if vepline[0].startswith('#'): # header
		for i in range(len(vepline)):
			index[vepline[i]] = i 
		continue
	if '_dupl' in vepline[index['Feature']]: # see https://www.biostars.org/p/312593/ , bug duplicate
		continue
	t = vepline[index['Feature']].split('.')[0]
	v = int(vepline[index['Feature']].split('.')[-1])
	if t in vep_transcript_version.keys():
		if v > vep_transcript_version[t]:
			vep_transcript_version[t] = v
vep_transcript_to_use = ['%s.%s' % (transcript,vep_transcript_version[transcript]) for transcript in vep_transcript_version.keys()]

vep_file.seek(0)
for vepline in vep_reader:
	if vepline[0].startswith('##'):
		continue
	if vepline[0].startswith('#'): # header
		continue
	location = vepline[index['Location']]
	chrom = 'chr%s' % location.split(':')[0]
	ref = vepline[index['GIVEN_REF']]
	alt = vepline[index['Allele']]
	if '-' in location :
		start = location.split(':')[-1].split('-')[0]
		stop = location.split(':')[-1].split('-')[-1]
	else:
		start = location.split(':')[-1]
		stop = start
	variant = '%s:%s-%s:%s>%s' % (chrom,start,stop,ref,alt)
	#print "-%s" % variant
	
	vep_transcript = vepline[index['Feature']]
	if not vep_transcript.startswith('NM_'): # Exclude NR, XM ...
		continue
	
	if variant not in annotation.keys():
		#print "\t-VEP Warning : VEP variant line does not match any ANNOVAR variant line"
		continue
	
	db_cur.execute("SELECT gene FROM Variant WHERE variantID='%s'" % variant)
	gene = db_cur.fetchone()['gene']
	exact_transcrit_with_version = '%s.%s' % (gene_data[gene]['transcript'],gene_data[gene]['transcriptVersion'])

	if vep_transcript in vep_transcript_to_use :
		if vep_transcript != exact_transcrit_with_version:
			annotation[variant]['annoWarning'].append("VEP transcript version : %s" % vep_transcript)
			#print "\t-VEP transcript version is %s" % vep_transcript
		if vepline[index['EXON']] != '-':
			annotation[variant]['vep_exon'] = representsInt(vepline[index['EXON']].split('/')[0])
		if vepline[index['INTRON']] != '-':
			annotation[variant]['vep_intron'] = representsInt(vepline[index['INTRON']].split('/')[0])
		if vepline[index['HGVSc']] != '-':
			annotation[variant]['vep_transcriptDescription'] = vepline[index['HGVSc']].split(':')[-1]
			if annotation[variant]['vep_transcriptDescription'] != annotation[variant]['transcriptDescription']:
				annotation[variant]['annoWarning'].append("VEP : %s" % annotation[variant]['vep_transcriptDescription'])
				#print "\t-VEP transcriptDescription is %s" % annotation[variant]['vep_transcriptDescription']
		if vepline[index['Consequence']] != '-':
			annotation[variant]['vep_consequence'] = vepline[index['Consequence']]
		if vepline[index['IMPACT']] != '-':
			annotation[variant]['vep_impact'] = vepline[index['IMPACT']]
		if vepline[index['PUBMED']] != '-':
			annotation[variant]['pubmed'] = vepline[index['PUBMED']]#.replace('-','.')
		vep_dbSNP = []
		vep_COSMIC = []
		if vepline[index['Existing_variation']] != '-':
			existing_variation = vepline[index['Existing_variation']].split(',')
			for item in existing_variation:
				if item.startswith('rs'):
					vep_dbSNP.append(item)
				elif item.startswith('COSM'):
					vep_COSMIC.append(item)
			if vep_dbSNP:
				annotation[variant]['vep_dbsnp'] = ','.join(vep_dbSNP)
			if vep_COSMIC:
				annotation[variant]['vep_cosmic'] = ','.join(vep_COSMIC)
		if vepline[index['CLIN_SIG']] != '-':
			annotation[variant]['vep_clinsig'] = vepline[index['CLIN_SIG']]
		if vepline[index['EA_AF']] != '-' and vepline[index['AA_AF']] != '-':
			try:
				annotation[variant]['vep_esp'] = str((float(vepline[index['EA_AF']]) + float(vepline[index['AA_AF']]))/2)
			except:
				annotation[variant]['vep_esp'] = None
		if vepline[index['AF']] != '-':
			annotation[variant]['vep_1000G_all'] = representsFloat(vepline[index['AF']])
		if vepline[index['EUR_AF']] != '-':
			annotation[variant]['vep_1000G_eur'] = representsFloat(vepline[index['EUR_AF']])
		if vepline[index['SIFT']] != '-':
			annotation[variant]['vep_sift'] = vepline[index['SIFT']]
		if vepline[index['PolyPhen']] != '-':
			annotation[variant]['vep_polyphen_hvar'] = vepline[index['PolyPhen']]

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
		if annotation[variant]['exon'] != None and annotation[variant]['vep_exon'] != None:
			if annotation[variant]['exon'] != annotation[variant]['vep_exon']:
				diff['exon'].append('+%s'%annotation[variant]['vep_exon'])
				diff['exon'].append('-%s'%annotation[variant]['exon'])
		elif annotation[variant]['exon'] == None and annotation[variant]['vep_exon'] != None:
			diff['exon'].append('+%s'%annotation[variant]['vep_exon'])
		elif annotation[variant]['exon'] != None and annotation[variant]['vep_exon'] == None:
			diff['exon'].append('-%s'%annotation[variant]['exon'])
				
		## CLINVAR DIFF (plusieurs valeurs possibles)
		annovar_clinsig = []
		vep_clinsig = []
		if annotation[variant]['clinsig'] != None:
			annovar_clinsig = annotation[variant]['clinsig'].replace('/',',').split(',')
		if annotation[variant]['vep_clinsig'] != None:
			vep_clinsig = annotation[variant]['vep_clinsig'].split(',')
		if 'not_provided' in vep_clinsig:
			vep_clinsig.remove('not_provided')
		for cs in vep_clinsig:
			if cs not in annovar_clinsig: # and cs != '-'
				diff['clinsig'].append('+%s'%cs)
		for cs in annovar_clinsig:
			if cs not in vep_clinsig: #  and cs != '.'
				diff['clinsig'].append('-%s'%cs)
				
		if annotation[variant]['clinsig'] == None and annotation[variant]['vep_clinsig'] != None: ## Si Annovar == None, utilisation VEP par defaut
			if 'benign' in annotation[variant]['vep_clinsig'] and 'pathogenic' in annotation[variant]['vep_clinsig']:
				annotation[variant]['clinsig'] = 'Conflicting_interpretations_of_pathogenicity'
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
				annotation[variant]['clinsig'] = ','.join(vep_clinsig)
				
		## COSMIC DIFF (plusieurs valeurs possibles)
		annovar_COSMIC = []
		if annotation[variant]['cosmic'] != None:
			annovar_COSMIC = annotation[variant]['cosmic'].split(',')
		for cosm in vep_COSMIC:
			if cosm not in annovar_COSMIC: #  and cosm != '.'
				diff['cosmic'].append('+%s'%cosm)
		for cosm in annovar_COSMIC:
			if cosm not in vep_COSMIC: #  and cosm != '.'
				diff['cosmic'].append('-%s'%cosm)
				
		## DBSNP DIFF (une valeur possible)
		if annotation[variant]['dbsnp'] != None and annotation[variant]['vep_dbsnp'] != None:
			if annotation[variant]['dbsnp'] != annotation[variant]['vep_dbsnp']:
				diff['dbsnp'].append('+%s'%annotation[variant]['vep_dbsnp'])
				diff['dbsnp'].append('-%s'%annotation[variant]['dbsnp'])
		elif annotation[variant]['dbsnp'] == None and annotation[variant]['vep_dbsnp'] != None:
			diff['dbsnp'].append('+%s'%annotation[variant]['vep_dbsnp'])
		elif annotation[variant]['dbsnp'] != None and annotation[variant]['vep_dbsnp'] == None:
			diff['dbsnp'].append('-%s'%annotation[variant]['dbsnp'])
			
		## 1000G_ALL DIFF (une valeur possible) # 2% de tolerance
		if annotation[variant]['1000G_all'] != None and annotation[variant]['vep_1000G_all'] != None:
			if not ((round(float(annotation[variant]['1000G_all']),2)-0.02) <= round(float(annotation[variant]['vep_1000G_all']),2) <= (round(float(annotation[variant]['1000G_all']),2)+0.02)):
				diff['1000G_ALL'].append('+%s'%annotation[variant]['vep_1000G_all'])
				diff['1000G_ALL'].append('-%s'%annotation[variant]['1000G_all'])
		elif annotation[variant]['1000G_all'] == None and annotation[variant]['vep_1000G_all'] != None:
			diff['1000G_ALL'].append('+%s'%annotation[variant]['vep_1000G_all'])
		elif annotation[variant]['1000G_all'] != None and annotation[variant]['vep_1000G_all'] == None:
			diff['1000G_ALL'].append('-%s'%annotation[variant]['1000G_all'])
			
		## 1000G_EUR DIFF (une valeur possible) # 2% de tolerance
		if annotation[variant]['1000G_eur'] != None and annotation[variant]['vep_1000G_eur'] != None:
			if not ((round(float(annotation[variant]['1000G_eur']),2)-0.02) <= round(float(annotation[variant]['vep_1000G_eur']),2) <= (round(float(annotation[variant]['1000G_eur']),2)+0.02)):
				diff['1000G_EUR'].append('+%s'%annotation[variant]['vep_1000G_eur'])
				diff['1000G_EUR'].append('-%s'%annotation[variant]['1000G_eur'])
		elif annotation[variant]['1000G_eur'] == None and annotation[variant]['vep_1000G_eur'] != None:
			diff['1000G_EUR'].append('+%s'%annotation[variant]['vep_1000G_eur'])
		elif annotation[variant]['1000G_eur'] != None and annotation[variant]['vep_1000G_eur'] == None:
			diff['1000G_EUR'].append('-%s'%annotation[variant]['1000G_eur'])
			
		## ESP DIFF (une valeur possible) # 2% de tolerance
		if annotation[variant]['esp'] != None and annotation[variant]['vep_esp'] != None:	
			if not ((round(float(annotation[variant]['esp']),2)-0.02) <= round(float(annotation[variant]['vep_esp']),2) <= (round(float(annotation[variant]['esp']),2)+0.02)):
				diff['esp'].append('+%s'%annotation[variant]['vep_esp'])
				diff['esp'].append('-%s'%annotation[variant]['esp'])
		elif annotation[variant]['esp'] == None and annotation[variant]['vep_esp'] != None:
			diff['esp'].append('+%s'%annotation[variant]['vep_esp'])
		elif annotation[variant]['esp'] != None and annotation[variant]['vep_esp'] == None:
			diff['esp'].append('-%s'%annotation[variant]['esp'])
			
		## SIFT DIFF (une valeur possible)	
		if annotation[variant]['sift'] != None and annotation[variant]['vep_sift'] != None:	
			if str(annotation[variant]['sift']).lower() != str(annotation[variant]['vep_sift']).replace('_',' ').replace('-','.'):
				diff['sift'].append('+%s'%annotation[variant]['vep_sift'])
				diff['sift'].append('-%s'%annotation[variant]['sift'])
		elif annotation[variant]['sift'] == None and annotation[variant]['vep_sift'] != None:
			diff['sift'].append('+%s'%annotation[variant]['vep_sift'])
		elif annotation[variant]['sift'] != None and annotation[variant]['vep_sift'] == None:
			diff['sift'].append('-%s'%annotation[variant]['sift'])
			
		## POLYPHEN DIFF (une valeur possible)
		if annotation[variant]['polyphen_hvar'] != None and annotation[variant]['vep_polyphen_hvar'] != None:	
			if str(annotation[variant]['polyphen_hvar']).lower() != str(annotation[variant]['vep_polyphen_hvar']).replace('_',' ').replace('-','.'):
				diff['polyphen_hvar'].append('+%s'%annotation[variant]['vep_polyphen_hvar'])
				diff['polyphen_hvar'].append('-%s'%annotation[variant]['polyphen_hvar'])
		elif annotation[variant]['polyphen_hvar'] == None and annotation[variant]['vep_polyphen_hvar'] != None:
			diff['polyphen_hvar'].append('+%s'%annotation[variant]['vep_polyphen_hvar'])
		elif annotation[variant]['polyphen_hvar'] != None and annotation[variant]['vep_polyphen_hvar'] == None:
			diff['polyphen_hvar'].append('-%s'%annotation[variant]['polyphen_hvar'])
			
		#### GENERATE DIFF STRING ####
		vep_diff = []
		for item in diff.keys():
			if diff[item]:
				d = ','.join(diff[item])
				vep_diff.append('%s:%s'%(item,d))
		if vep_diff:
			 annotation[variant]['vep_diff'] = ';'.join(vep_diff)#annotation[variant]['vep_diff'][1:]
		
########################################
## COMPARING ANNOVAR/VEP CONSEQUENCES ##
########################################

print "- [%s] comparing Annovar/VEP consequences ..." % (time.strftime("%H:%M:%S"))
for variant in annotation.keys():
	if annotation[variant]['consequence'] == None:
		if annotation[variant]['vep_consequence'] != None:
			vep_consequences = annotation[variant]['vep_consequence'].split(',') # sometimes several consequences
			region_consequence_list = []
			for consequence in vep_consequences:
				if consequence in vep2annovar:
					region_consequence_list.append((vep2annovar[consequence]['region'],vep2annovar[consequence]['type']))
			final_consequence = compare_consequence(region_consequence_list)
			annotation[variant]['region'] = final_consequence[0]
			annotation[variant]['consequence'] = final_consequence[1]
			annotation[variant]['dbsnp'] = annotation[variant]['vep_dbsnp']
			annotation[variant]['cosmic'] = annotation[variant]['vep_cosmic']
			annotation[variant]['dbsnp'] = annotation[variant]['vep_dbsnp']
			annotation[variant]['clinsig'] = annotation[variant]['vep_clinsig']
			annotation[variant]['esp'] = annotation[variant]['vep_esp']
			annotation[variant]['1000G_all'] = annotation[variant]['vep_1000G_all']
			annotation[variant]['sift'] = annotation[variant]['vep_sift']
			annotation[variant]['polyphen_hvar'] = annotation[variant]['vep_polyphen_hvar']
			annotation[variant]['annoWarning'].append("VEP annotation was used")
			#annotation[variant]['commentaire'].append("VEP annotation")
			#print "- compare %s : VEP consequence and annotations were used to annotate this variant" % variant
	elif annotation[variant]['vep_consequence'] != None:
		vep_consequences = annotation[variant]['vep_consequence'].split(',') # sometimes several consequences
		region_consequence_list = [(annotation[variant]['region'],annotation[variant]['consequence'])]
		for consequence in vep_consequences:
			if consequence in vep2annovar:
				region_consequence_list.append((vep2annovar[consequence]['region'],vep2annovar[consequence]['type']))
		final_consequence = compare_consequence(region_consequence_list)
		if final_consequence != (annotation[variant]['region'],annotation[variant]['consequence']):
			annotation[variant]['annoWarning'].append("Annovar/VEP discordance for consequence")
			#print "- compare %s : Annovar/VEP discordance for consequence" % variant
	# TODO : modifier vep2annovar et enlever les '.', en attendant:
	if annotation[variant]['consequence'] == '.':
		annotation[variant]['consequence'] = None


#############################
# ADDING CUSTOM ANNOTATIONS #
#############################

# KNOWN FALSE-POSITIVES, DRUG SENSITIVITY, LAM HOTSPOT HIGHLIGHTS, TP53 UMD
print "- [%s] adding custom annotations ..." % (time.strftime("%H:%M:%S"))

custom_fp = []
with open(custom_false_positives,'r') as fp:
	for row in fp :
		fp_row = row.replace('\n','').split(';')
		custom_fp.append(fp_row)

sensi = []
sensitivity_xls = xlrd.open_workbook(custom_drug_sensitivity)
sensitivity_sheet = sensitivity_xls.sheet_by_index(0)
for row in range(sensitivity_sheet.nrows-1):
	sensi.append((str(sensitivity_sheet.cell(row+1,5).value),str(sensitivity_sheet.cell(row+1,7).value),str(sensitivity_sheet.cell(row+1,8).value),str(sensitivity_sheet.cell(row+1,9).value))) # gene,c,p,sensi

custom_hl = []
highlight_xls = xlrd.open_workbook(custom_highlight)
highlight_sheet = highlight_xls.sheet_by_index(0)
for row in range(highlight_sheet.nrows-1):
	custom_hl.append((highlight_sheet.cell(row+1,5).value,highlight_sheet.cell(row+1,7).value)) # ,highlight_sheet.cell(row+1,8).value
	
custom_tu = []
with open(custom_tp53_umd,'r') as tp53_umd:
	tp53_umd.next()
	for row in tp53_umd:
		tp53_row = row.replace('\n','')
		tp53_row = tp53_row.split('\t')
		custom_tu.append((tp53_row[0],tp53_row[1],tp53_row[2],tp53_row[3]))
		
custom_lc = []
with open(custom_lymphome_controls,'r') as lc:
	for row in lc :
		lc_row = row.replace('\n','').split('\t')
		nb = len(lc_row[2].split(','))
		string = "Found in normal dna (x%s)" % nb # "Found in %s controls (%s)" % (nb,lc_row[2])
		lc_row[2] = string
		custom_lc.append(lc_row)

for variant in annotation.keys():
	# KNOWN FALSE-POSITIVES
	for fp in custom_fp:
		if (fp[0] == annotation[variant]['transcript']) and ((fp[1] == annotation[variant]['transcriptDescription']) or (fp[1] == annotation[variant]['annovar_transcriptDescription'])):
			annotation[variant]['commentaire'].append(fp[2])
			break
	# DRUG SENSITIVITY
	for s in sensi:
		if s[2] not in ['','p.(=)','p.?','p.=']:
			if ((annotation[variant]['gene'],annotation[variant]['proteinDescription']) == (s[0],s[2])) or ((annotation[variant]['gene'],annotation[variant]['annovar_proteinDescription']) == (s[0],s[2])) :
				annotation[variant]['actionability'] = s[3]
				break
		if s[1] != '':
			if ((annotation[variant]['gene'],annotation[variant]['transcriptDescription']) == (s[0],s[1])) or ((annotation[variant]['gene'],annotation[variant]['annovar_transcriptDescription']) == (s[0],s[1])):
				annotation[variant]['actionability'] = s[3]
				break
	## cas particulier : splicing MET
	if (annotation[variant]['gene'] == 'MET') and (annotation[variant]['region'] == 'intronic') and ((116411873 <= annotation[variant]['genomicStart'] <= 116411902) or (116412044 <= annotation[variant]['genomicStart'] <= 116412087)):
		annotation[variant]['exon'] = 14
		annotation[variant]['commentaire'].append('MET intron 13/14')
		if annotation[variant]['actionability'] == None :
			annotation[variant]['actionability'] = 'intron 13/14'
	# LAM HOTSPOT HIGHLIGHTS
	for hl in custom_hl:
		if ((annotation[variant]['gene'],annotation[variant]['transcriptDescription']) == hl) or ((annotation[variant]['gene'],annotation[variant]['annovar_transcriptDescription']) == hl): # or ((annotation[variant]['gene'],annotation[variant]['transcriptDescription'],annotation[variant]['proteinDescription']) == hl) # necessite hgvs avant
			annotation[variant]['highlight'] = 1
			break
	#TP53 UMD
	for tu in custom_tu:
		if tu[1] not in ['','p.(=)','p.?','p.=']:		
			if ((annotation[variant]['gene'],annotation[variant]['proteinDescription']) == ('TP53',tu[1])) or ((annotation[variant]['gene'],annotation[variant]['annovar_proteinDescription']) == ('TP53',tu[1])):
				annotation[variant]['pathoUMD'] = tu[2]
				annotation[variant]['commentUMD'] = tu[3]
				break
		if tu[0] not in ['','c.(=)','c.?','c.=']:		
			if ((annotation[variant]['gene'],annotation[variant]['transcriptDescription']) == ('TP53',tu[0])) or ((annotation[variant]['gene'],annotation[variant]['annovar_transcriptDescription']) == ('TP53',tu[0])):
				annotation[variant]['pathoUMD'] = tu[2]
				annotation[variant]['commentUMD'] = tu[3]
				break
	#LYMPHOME CONTROLS
	for lc in custom_lc:
		if (lc[0] == annotation[variant]['transcript']) and ((lc[1] == annotation[variant]['transcriptDescription']) or (lc[1] == annotation[variant]['annovar_transcriptDescription'])):
			annotation[variant]['commentaire'].append(lc[2])
			break
	
#################
# WRITE RESULTS #
#################	

print "- [%s] writting results to VariantBase ..." % (time.strftime("%H:%M:%S"))

lastUpdate = '%s' % (time.strftime('%d/%m/%Y'))
for variant in annotation.keys():
	if annotation[variant]['commentaire']:
		commentaire = '. '.join(annotation[variant]['commentaire'])
	else:
		commentaire = None
	if annotation[variant]['annoWarning']:
		annoWarning = '. '.join(annotation[variant]['annoWarning'])
	else:
		annoWarning = None

	try:
		db_cur.execute("""UPDATE Variant SET exon=?,intron=?,region=?,consequence=?,actionability=?,intervar=?,clinvar=?,cosmic=?,dbsnp=?,gnomad=?,milleGall=?,milleGeur=?,nci60=?,esp=?,exac=?,sift=?,polyphen2=?,provean=?,pubmed=?,vep_consequence=?,vep_impact=?,vep_diff=?,commentaire=?,highlight=?,pathoUMD=?,commentUMD=?,annovarTranscriptDescription=?,annovarProteinDescription=?,annoWarning=?,lastUpdate=? WHERE variantID=?""", (annotation[variant]['exon'],annotation[variant]['vep_intron'],annotation[variant]['region'],annotation[variant]['consequence'],annotation[variant]['actionability'],annotation[variant]['intervar'],annotation[variant]['clinsig'],annotation[variant]['cosmic'],annotation[variant]['dbsnp'],annotation[variant]['gnomad'],annotation[variant]['1000G_all'],annotation[variant]['1000G_eur'],annotation[variant]['nci60'],annotation[variant]['esp'],annotation[variant]['exac'],annotation[variant]['sift'],annotation[variant]['polyphen_hvar'],annotation[variant]['provean'],annotation[variant]['pubmed'],annotation[variant]['vep_consequence'],annotation[variant]['vep_impact'],annotation[variant]['vep_diff'],commentaire,annotation[variant]['highlight'],annotation[variant]['pathoUMD'],annotation[variant]['commentUMD'],annotation[variant]['annovar_transcriptDescription'],annotation[variant]['annovar_proteinDescription'],annoWarning,lastUpdate,variant))
	except Exception as e:
		print "Error : %s" % e

db_con.commit()
db_con.close()
print "- [%s] ... done." % (time.strftime("%H:%M:%S"))
