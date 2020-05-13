#!/usr/bin/python
import sqlite3
import zipfile
import json
import uuid
import sys
import os
from datetime import date

# THIS SCRIPT FILL AN EMPTY VariantBase.db
# IT USE A RUN LIST AND VARIANT LIST TO CREATE AND FILL ALL TABLES : RUN, SAMPLE, ANALYSIS, VARIANT, VARIANTMETRICS
# (USE insert_all_panel_variantbase.py FIRST FOR PANEL AND TARGETEDREGIONS)

# USAGE : python insert_all_run_all_samples_variantbase.py

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))
	
db_path = global_param['VariantBase']

month2num = {'Jan':1,'Feb':2,'Mar':3,'Apr':4,'May':5,'Jun':6,'Jul':7,'Aug':8,'Sep':9,'Oct':10,'Nov':11,'Dec':12}
run_list_path = '%s/variantBase/runList/runList_ALL.fullpath.txt' % pipeline_folder
variant_list_path = '%s/variantBase/variantList/variantList_ALL.json' % pipeline_folder

with open(variant_list_path, 'r') as g:
	print "- loading variantList_ALL.json..."
	variantlist = json.load(g)

#   ___                     __      __               ___       __        ___ 
#  |__  | |    |    | |\ | / _`    |__) |  | |\ |     |   /\  |__) |    |__  
#  |    | |___ |___ | | \| \__>    |  \ \__/ | \|     |  /~~\ |__) |___ |___ 	

run_list = []
rl = open(run_list_path,'r')
for run in rl:
	run_folder = run.replace('\n','')
	run_list.append(run_folder)
		
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

db_cur.execute("SELECT * FROM Gene")
db_genes = db_cur.fetchall()
gene_data = {}
for db_gene in db_genes:
	gene_data[db_gene['geneID']] = {
		'chromosome': db_gene['chromosome'],
		'transcriptionStart': db_gene['transcriptionStart'],
		'transcriptionStop': db_gene['transcriptionStop']
	}

for run_folder in run_list:
	# RUN TABLE
	if run_folder.endswith('/'):
		run_name = run_folder.split('/')[-2]
	else:
		run_name = run_folder.split('/')[-1]
		
	db_cur.execute("SELECT runID FROM Run WHERE runID='%s'" % run_name)
	if db_cur.fetchone():
		print "* Run is already in database (%s)" % run_name
		continue

	if 'S5' in run_name: # WARNING : make sure run name all have PGM or S5 in their name
		platform = 'S5'
	if 'PGM' in run_name:
		platform = 'PGM'

	archive = zipfile.ZipFile(run_folder+'/pgm_logs.zip', 'r')

	with archive.open('explog_final.txt') as explog_final:
		for line in explog_final:
			if line.startswith('Start Time:'):
				start_time = line.split('Start Time:')[-1].strip()
				if platform == 'PGM': # start time look like :  Tue Jan 24 08:54:47 2017
					start_time = start_time.replace('  ',' ')
					if start_time.startswith(' '):
						start_time = start_time[1:]
					start_time = start_time.split(' ')
					day = int(start_time[2])
					month = month2num[start_time[1]]
					year = int(start_time[4])
					run_date = date(year,month,day)
				elif platform == 'S5': # start time look like :  01/18/2017 10:23:51
					start_time = start_time.replace(' ','/').split('/')
					day = int(start_time[1])
					month = int(start_time[0])
					year = int(start_time[2])
					run_date = date(year,month,day)
				break
			else:
				continue
				
	print "#######################################"
	print "RUN : %s" % run_name # PK
	print "PLATFORM : %s" % platform
	print "DATE : %s" % run_date
	print "FOLDER : %s" % run_folder
	print "#######################################"

	try:
		db_cur.execute("INSERT INTO Run (runID, platform, runPath, runDate) VALUES ('%s', '%s', '%s', '%s')" % (run_name,platform,run_folder,run_date))
	except Exception as e:
		print "*WARNING* (RUN table)** %s" % e

	#   ___                     __      __              __        ___    ___       __        ___ 
	#  |__  | |    |    | |\ | / _`    /__`  /\   |\/| |__) |    |__      |   /\  |__) |    |__  
	#  |    | |___ |___ | | \| \__>    .__/ /~~\  |  | |    |___ |___     |  /~~\ |__) |___ |___ 
                                                                                         
	with open(run_folder+'/barcodes.json', 'r') as g:
		barcodes_json = json.load(g)
	for barcode in barcodes_json:
		sample = barcodes_json[barcode]['sample']
		dna_number = barcodes_json[barcode]['sample_id']
		panel = barcodes_json[barcode]['target_region_filepath'].split('/unmerged/detail/')[-1]
		#tumorCellularity = "?" # voir bilan temoin (optionnel)
		#gender = "?" # voir bilan temoin (optionnel)
		iscontrol = 0
		if dna_number.startswith('CONTROL-'):
			iscontrol = 1
		name = sample.replace(dna_number,'')
		while name.endswith('-') or name.endswith('_'):
			name = name[:-1]
		print "- %s :" % sample
		print "\t- dna_number -> %s" % dna_number # PK
		print "\t- name -> %s" % name
		
		db_cur.execute("SELECT sampleID FROM Sample WHERE sampleID='%s'" % dna_number)
		if db_cur.fetchone() is None:
			try:
				db_cur.execute("INSERT INTO Sample (sampleID, sampleName, isControl) VALUES ('%s', '%s', '%s')" % (dna_number,name,iscontrol))
			except Exception as e:
				print "\t*WARNING* (SAMPLE table)** %s" % e
		else:
			print "\t* (sample is already in database)"
			
		#   ___                     __                              __     __     ___       __        ___ 
		#  |__  | |    |    | |\ | / _`     /\  |\ |  /\  |    \ / /__` | /__`     |   /\  |__) |    |__  
		#  |    | |___ |___ | | \| \__>    /~~\ | \| /~~\ |___  |  .__/ | .__/     |  /~~\ |__) |___ |___ 
																							
		random_uuid = uuid.uuid1()
		analysis_id = 'A-'+random_uuid.hex[:8]
		bam_path = '%s/%s/%s_%s.bam' % (run_folder,sample,sample,barcode)
		db_cur.execute("SELECT panelID FROM Panel WHERE panelID='%s'"%panel)
		if db_cur.fetchone() is None:
			print "\t*WARNING* (PANEL %s not found)** " % panel
			print "\t** Panel is needed for foreign key constraint. This sample will not be included."
			continue
		if not os.path.exists(bam_path):
			print "\t*WARNING* (BAM PATH %s inexistant )** " % bam_path
			bam_path = ''
		try:
			db_cur.execute("INSERT INTO Analysis (analysisID, sample, barcode, run, panel, bamPath, analysisDate) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s')" % (analysis_id, dna_number, barcode, run_name, panel, bam_path, run_date))
		except Exception as e:
			print "\t*WARNING* (ANALYSIS table)** %s"%e
			
		#   ___                     __                __              ___    ___       __        ___ 
		#  |__  | |    |    | |\ | / _`    \  /  /\  |__) |  /\  |\ |  |      |   /\  |__) |    |__  
		#  |    | |___ |___ | | \| \__>     \/  /~~\ |  \ | /~~\ | \|  |      |  /~~\ |__) |___ |___ 
		
		genome_build = 'hg19'
		vc_tool = 'TVC'		
		if run_name in variantlist:
			if sample in variantlist[run_name]:
				print "\t- %s variants" % len(variantlist[run_name][sample])
				for variant in variantlist[run_name][sample]:
					chromosome = variant[0]
					start = variant[1]
					stop = variant[2]
					ref = variant[3]
					alt = variant[4]
					variant_id = '%s:%s-%s:%s>%s' % (chromosome,start,stop,ref,alt)
					varcov = variant[5]
					poscov = variant[6]
					
					#Genomic Description
					if ref == '-':
						genomicDescription = '%s:g.%s_%sins%s' % (chromosome,start,stop,alt)
					elif alt == '-':
						if len(ref) > 1:
							genomicDescription = '%s:g.%s_%sdel%s' % (chromosome,start,stop,ref)
						else:
							genomicDescription = '%s:g.%sdel%s' % (chromosome,start,ref)
					elif len(ref) > 1 or len(alt) > 1:
						if len(ref) > 1:
							genomicDescription = '%s:g.%s_%sdelins%s' % (chromosome,start,stop,alt)
						else:
							genomicDescription = '%s:g.%sdelins%s' % (chromosome,start,alt)
					else:
						genomicDescription = '%s:g.%s%s>%s' % (chromosome,start,ref,alt)
					
					gene = ''
					for g in gene_data.keys():
						if chromosome == gene_data[g]['chromosome'] and ((gene_data[g]['transcriptionStart']-5000)<int(start)<(gene_data[g]['transcriptionStop']+5000)):
							gene = g

					if ref == '-': 
						variant_type = 'INS'
					elif alt == '-': 
						variant_type = 'DEL'
					elif len(ref) > 1 or len(alt) > 1: 
						variant_type = 'DELINS'
					else:
						variant_type = 'SNV'
					
					## EXON ?
					#for index in range(len(exonStarts)):
						#if (int(exonStarts[index]) < int(ampl_start) < int(exonEnds[index])) or (int(exonStarts[index]) < int(ampl_stop) < int(exonEnds[index])):
							#if strand == 'forward':
								#exon = index + 1
							#elif strand == 'reverse':
								#exon = len(exonStarts)-index
							#break
					## INTRON ?
					#for index in range(len(exonStarts)-1):
						#if (int(exonEnds[index]) < int(ampl_start) < int(exonStarts[index+1])) or (int(exonEnds[index]) < int(ampl_stop) < int(exonStarts[index+1])):
							#if strand == 'forward':
								#intron = index + 1
							#elif strand == 'reverse':
								#intron = (len(exonStarts)-index)-1
							#break
					
					db_cur.execute("SELECT variantID FROM Variant WHERE variantID='%s'" % variant_id)
					if db_cur.fetchone() is None:
						try:
							#print "\t\t * new variant : %s)" % variant
							db_cur.execute("INSERT INTO Variant (variantID, genomeBuild, chromosome, genomicStart, genomicStop, referenceAllele, alternativeAllele, variantType, gene, genomicDescription) VALUES ('%s','%s','%s',%s, %s,'%s','%s','%s','%s','%s')" % (variant_id, genome_build, chromosome, start, stop, ref, alt,variant_type,gene,genomicDescription))
						except Exception as e:
							print "\t*WARNING* (VARIANT table)** %s"%e
							
					#   ___                     __                __              ___        ___ ___  __     __   __     ___       __        ___ 
					#  |__  | |    |    | |\ | / _`    \  /  /\  |__) |  /\  |\ |  |   |\/| |__   |  |__) | /  ` /__`     |   /\  |__) |    |__  
					#  |    | |___ |___ | | \| \__>     \/  /~~\ |  \ | /~~\ | \|  |   |  | |___  |  |  \ | \__, .__/     |  /~~\ |__) |___ |___ 
																															  
					random_uuid = uuid.uuid1()
					variantmetrics_id = 'M-'+random_uuid.hex[:8]
					
					try:
						db_cur.execute("INSERT INTO VariantMetrics (variantMetricsID, variant, analysis, positionReadDepth, variantReadDepth, variantCallingTool, call) VALUES ('%s', '%s', '%s', %s, %s, '%s', 'de novo')" % (variantmetrics_id, variant_id, analysis_id, poscov, varcov, vc_tool))
					except Exception as e:
						print "\t*WARNING* (VARIANTMETRICS table)** %s"%e
			else:
				print "\t*!WARNING!* (SAMPLE MISSING RESULTS)** "
		else:
			print "\t*!WARNING!* (RUN MISSING RESULTS)** "


##########################################################################################################

#   __   ___        __               __      __         ___ ___                   __              ___  __  
#  |__) |__   |\/| /  \ \  / | |\ | / _`    /__` |__| |  |   |  \ /    \  /  /\  |__) |  /\  |\ |  |  /__` 
#  |  \ |___  |  | \__/  \/  | | \| \__>    .__/ |  | |  |   |   |      \/  /~~\ |  \ | /~~\ | \|  |  .__/ 
                                                                                                        
# bad variant to correct
db_cur.execute("SELECT * FROM VariantMetrics WHERE variant='chr7:55242469-55242486:TAAGAGAAGCAACATCTC>-'")
varmetrics = db_cur.fetchall()
for varmetric in varmetrics:
	VariantMetricsID = varmetric['variantMetricsID']
	db_cur.execute("UPDATE VariantMetrics SET variant='chr7:55242469-55242486:TTAAGAGAAGCAACATCT>-' WHERE VariantMetricsID = '%s'" % VariantMetricsID)
	
shit_variants = [
'chr7:55242469-55242486:TAAGAGAAGCAACATCTC>-',
'chr9:133729624-133729627:GGTG>AGGTGA',
'chr9:133738339-133738419:AAGCTGGGCGGGGGCCAGTACGGGGAGGTGTACGAGGGCGTGTGGAAGAAATACAGCCTGACGGTGGCCGTGAAGACCTTG>-',
'chr9:133748422-133748493:CAGAGATCTTGCTGCCCGAAACTGCCTGGTAGGGGAGAACCACTTGGTGAAGGTAGCTGATTTTGGCCTGAG>-',
'chr11:534331-534440:AGGGGCCTGCGGCCCGGGGTCCTCCTACAGGGTCTCCTGCCCCACCTGCCAAGGAGGGCCCTGCTCAGCCAGGCCCAGGCCCAGCCCCAGGCCCCACAGGGCAGCTGCTG>-',
'chr7:116411842-116411914:AAGTCTCCTGGGGCCCATGATAGCCGTCTTTAACAAGCTCTTTCTTTCTCTCTGTTTTAAGATCTGGGCAGTG>-',
'chr1:2489255-2489358:TGCCCCAAGTGCAGTCCAGGTAGGTGCAGCCCTTTGGCGGGCCAGCTCTGTGGGCCGAGGGCAGACACTCTTGCCCCCTTCTGCCCCAGACACCCCTGTGTTCT>-',
'chr10:31816312-31816382:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAA>AGTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACAAC',
'chr10:31816312-31816382:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAA>AGTAAAAACTAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACAAC',
'chr10:31816312-31816381:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACA>AGTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACAAC',
'chr17:7753023-7753135:TGCTGACGTCGTGCGCGCCAGCAGGTGAGTCGGCTGCCTGCTTGCTTGTCCCGGAGACAGGCTCCCTTCCCCCATCACCCTGATGTTTCTGTCTTCATTCCACTGTCCTCCCG>CGCCGATGTCATGCGTGCCAGCAAGTGAGTTCGGAGGCTGCCGCACTTGTCCTGGGGCTAGGCTTGTGCTCGTTACCCCATGTTCCTGACCTCATTCTATGGTCTTCTCA',
'chr17:7751226-7751329:TCACACCCCTCCCACTCCCCCAACCCCAACCACCAGCAGTAGCAACAGCAACAGTGGCAGCCACAGCAGCAGCCCTGCTGGGCCTGTGTCCTTTCCCCCACCAC>GCATAACCCTCCCATTCCCCCAACCACCACCAGCAGCAGCAGCAGCAGCAACAGCCACAGCAGTAGTCCTACTGGGCCGGTGCCCTTTCCACCACCCT',
'chr10:31816312-31816381:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACA>AGTAAAAACTAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACAAC',
'chr10:31816312-31816382:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAA>AGTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAAC',
'chr10:31816312-31816381:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACA>AGTAAAAAACTAAAAAAATACAAAAATACAAAACACACACACACACACACACACACACACACACACACACAAC',
'chr10:31816312-31816382:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAA>AGTAAAAACTAAAAAATACAAAAATACAAAACACACACACACACACACACACACACACACACACAAC',
'chr10:31816319-31816382:CTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAA>ACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACAAC',
'chr20:39788406-39788500:TTAAGGGGGTAGAGGAGGTAGAGGATAGTTAGGGGAATGCCTGCTGGCTCCTGCCCAGTGGGAGGTATGTGCCCTCGGGGCAGCTATTGATACCT>ACATGGGCTGAGAGGAGGTGGGGTGTAGGTAAGGGGAGCCCAAGCACTTGCTCTTGTGGGAGCTGTGTGCCCTTGGTTCAGCCACTGCTC',
'chr10:31816312-31816381:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACA>AGTAAAAACTAAAAAAATACAAAAATACAAAACACACACACACACACACACACACACACACACACACACAAC',
'chr10:31816312-31816382:GTAAAAACTAAAAAAATACAAAATACAAAACACACACACACACACACACACACACACACACACACACACAA>AGTAAAAACTAAAAAAATACAAAAATACAAAACACACACACACACACACACACACACACACACACAAC',
'chr9:133738420-133738421:->AGGGGAGCTGCTGGTGAGGATTATTTTAGACTGTGAGTAATTGACCTGACAGACAGTGATGACTGCTTCATTAAGAGCCCACGACCACGTGCCAGAATAGTTCAGCATCCTCTGTTGCTACTGTACTTTGAGACATCGTTCTTCTTTGTGATGCAATACCTCTTTCTTGTCATGAGGGTCTCTTCCCTTAAATC',
'chr9:133748422-133748423:->AGTCTGAATAGGACAAGTGGAGATGCACTCTGATAGAAGTTTCTTGCAATCAAGGACTGGGCTCAGTGTTTGTCTCCGTGCCTGGCATAGAGTAGCAATAACTGCACTTACTATTCTCGACACTTCAATCAAAGCACTTAGCATATTAATTAATTACACTCACAGCAACCCTAAG'
]
for sv in shit_variants:
	db_cur.execute("SELECT * FROM VariantMetrics WHERE variant='%s'" % sv)
	varmetrics = db_cur.fetchall()
	for varmetric in varmetrics:
		db_cur.execute("DELETE FROM VariantMetrics WHERE VariantMetricsID='%s'" % varmetric['variantMetricsID'])
	db_cur.execute("DELETE FROM Variant WHERE variantID='%s'" % sv)
				
db_con.commit()
db_con.close()
