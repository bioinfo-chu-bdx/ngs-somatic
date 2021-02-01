#!/usr/bin/python
import sqlite3
import openpyxl
import pysam
import glob
import json
import csv
import os

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
db_con = sqlite3.connect('%s/variantBase/VariantBase.db' % pipeline_folder)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

 # __   ___       ___ ___  ___     ___       ___  __      ___              __      __   ___           ___  ___  __     ___  __              __                                            __     __                  __              ___        ___ ___  __     __   __                __      __             
# |  \ |__  |    |__   |  |__     |__  \  / |__  |__) \ /  |  |__| | |\ | / _`    |__) |__  |     /\   |  |__  |  \     |  /  \     /\     |__) |  | |\ |    .    /\  |\ |  /\  |    \ / /__` | /__`      \  /  /\  |__) |  /\  |\ |  |   |\/| |__   |  |__) | /  ` /__`     /\  |\ | |  \    |__) |  | |\ |  
# |__/ |___ |___ |___  |  |___    |___  \/  |___ |  \  |   |  |  | | | \| \__>    |  \ |___ |___ /~~\  |  |___ |__/     |  \__/    /~~\    |  \ \__/ | \|    .   /~~\ | \| /~~\ |___  |  .__/ | .__/ .     \/  /~~\ |  \ | /~~\ | \|  |   |  | |___  |  |  \ | \__, .__/    /~~\ | \| |__/    |  \ \__/ | \| .

# run = 'Auto_user_S5-0198-470-NGS518_519_35pM_Chef_SBT-colon-lung_v10_530_700'
# db_cur.execute("DELETE FROM VariantMetrics WHERE variantMetricsID IN (SELECT variantMetricsID FROM VariantMetrics INNER JOIN Analysis ON Analysis.analysisID = VariantMetrics.analysis INNER JOIN Run ON Run.runID = Analysis.run WHERE Run.runID='%s')" % run)
# db_cur.execute("DELETE FROM Analysis WHERE analysisID IN (SELECT analysisID FROM Analysis INNER JOIN Run ON Run.runID = Analysis.run WHERE Run.runID='%s')" % run)
# db_cur.execute("DELETE FROM Run WHERE Run.runID='%s'" % run)


        # __   ___  __  ___     ___      __   __   __  ___  ___  __      __   __      __       ___      
# | |\ | /__` |__  |__)  |     |__  \_/ |__) /  \ |__)  |  |__  |  \    |  \ |__)    |  \  /\   |   /\  
# | | \| .__/ |___ |  \  |     |___ / \ |    \__/ |  \  |  |___ |__/    |__/ |__)    |__/ /~~\  |  /~~\ 
                                                                                                      
	# z# # # # SELECT * FROM Analysis
	# # # # # INNER JOIN Run ON Run.runID = Analysis.run
	# # # # # INNER JOIN Sample ON Sample.sampleID = Analysis.sample
	# # # # # INNER JOIN Panel ON Panel.panelID = Analysis.panel
	# # # # # INNER JOIN VariantMetrics ON VariantMetrics.analysis = Analysis.analysisID
	# # # # # WHERE runID in ('Auto_user_S5-0198-507-NGS540_541_35pM_Chef_SBT-colon-lung_v10_530_738','Auto_user_S5-0198-508-Run224-TP53-ABL-MM_739',....

# added_runs = []
# added_samples = []
# added_analysis = []
# added_variants = []
# with open('/media/stuff/variantBase_export/export.csv','r') as export:
	# export_reader = csv.DictReader(export,delimiter=';')
	# for line in export_reader:
		# if line['runID'] not in added_runs:
			# db_cur.execute("INSERT INTO Run (runID,platform,system,runPath,runDate) VALUES ('%s','%s','%s','%s',%s)" % (line['runID'],line['platform'],line['system'],line['runPath'],int(line['runDate'])))
			# added_runs.append(line['runID'])
		# if line['sampleID'] not in added_samples:
			# try:
				# db_cur.execute("INSERT INTO Sample (sampleID,sampleName,gender,pathology,isControl) VALUES ('%s','%s','%s','%s',%s)" % (line['sampleID'],line['sampleName'],line['gender'],line['pathology'],int(line['isControl'])))
			# except Exception as e:
				# print str(e) + line['sampleID']
			# added_samples.append(line['sampleID'])
		# if line['analysisID'] not in added_analysis:
			# db_cur.execute("INSERT INTO Analysis (analysisID,sample,barcode,run,panel,bamPath,analysisDate) VALUES ('%s','%s','%s','%s','%s','%s',%s)" % (line['analysisID'],line['sample'],line['barcode'],line['run'],line['panel'],line['bamPath'],int(line['analysisDate'])))
			# added_analysis.append(line['analysisID'])
		# if line['variant'] not in added_variants:
			# db_cur.execute("INSERT OR IGNORE INTO Variant (variantID,chromosome,genomicStart,genomicStop,referenceAllele,alternativeAllele) VALUES ('%s','%s',%s,%s,'%s','%s')" % (line['variant'],line['chromosome'],int(line['genomicStart']),int(line['genomicStop']),line['referenceAllele'],line['alternativeAllele']))
			# added_variants.append(line['variant'])
		# db_cur.execute("INSERT INTO VariantMetrics (variantMetricsID,analysis,variant,positionReadDepth,variantReadDepth,variantCallingTool,call) VALUES ('%s','%s','%s',%s,%s,'%s','%s')" % (line['variantMetricsID'],line['analysis'],line['variant'],int(line['positionReadDepth']),int(line['variantReadDepth']),line['variantCallingTool'],line['call']))
		# # db_cur.execute("UPDATE VariantMetrics SET call='%s' WHERE variantMetricsID='%s'" % (line['call'],line['variantMetricsID']))

# db_con.commit()
# db_con.close()

 # __                  __   ___     __                     ___  __   __             ___ 
# /  ` |__|  /\  |\ | / _` |__     /  `  /\  |    |       |__  /  \ |__)  |\/|  /\   |  
# \__, |  | /~~\ | \| \__> |___    \__, /~~\ |___ |___    |    \__/ |  \  |  | /~~\  |  
                                                                                      

# db_cur.execute("SELECT * FROM VariantMetrics")
# db_trucs = db_cur.fetchall()
# for db_truc in db_trucs:
	# corrected_truc = False
	# truc_id = db_truc['variantMetricsID']
	# call = db_truc['call']
	# if db_truc['variantCallingTool'] == 'TVC':
		# corrected_truc = call.replace(' / ','/')
		# corrected_truc = corrected_truc.replace('de novo','tvc_de_novo')
		# corrected_truc = corrected_truc.replace('hotspot','tvc_hotspot')
	# if corrected_truc:
		# db_cur.execute("UPDATE VariantMetrics SET call='%s' WHERE variantMetricsID='%s'" % (corrected_truc,truc_id))

 # __                  __   ___     __       ___  ___     ___  __   __             ___ 
# /  ` |__|  /\  |\ | / _` |__     |  \  /\   |  |__     |__  /  \ |__)  |\/|  /\   |  
# \__, |  | /~~\ | \| \__> |___    |__/ /~~\  |  |___    |    \__/ |  \  |  | /~~\  |  
                                                                                     

# db_cur.execute("SELECT * FROM VariantAnnotation")
# db_trucs = db_cur.fetchall()
# for db_truc in db_trucs:
	# truc_id = db_truc['variantAnnotationID']
	# truc_date = db_truc['lastUpdate']
	# if '/' in truc_date:
		# truc_date = truc_date.split('/')
		# j = truc_date[0]
		# m = truc_date[1]
		# a = truc_date[2]
	# elif '-' in truc_date:
		# truc_date = truc_date.split('-')
		# j = truc_date[2]
		# m = truc_date[1]
		# a = truc_date[0]
	# corrected_truc_date = '%s%s%s' % (a,m,j)
	# db_cur.execute("UPDATE VariantAnnotation SET lastUpdate='%s' WHERE variantAnnotationID='%s'" % (corrected_truc_date,truc_id))

 # ___         __                __      __   __   __   __   ___  __  ___     __       ___                    ___          __                  __            __        ___  __  
# |__  | |\ | |  \     /\  |\ | |  \    /  ` /  \ |__) |__) |__  /  `  |     |__)  /\   |  |__|    |  | |__| |__  |\ |    |__) |  | |\ |    | /__`     |\/| /  \ \  / |__  |  \ 
# |    | | \| |__/    /~~\ | \| |__/    \__, \__/ |  \ |  \ |___ \__,  |     |    /~~\  |  |  |    |/\| |  | |___ | \|    |  \ \__/ | \|    | .__/     |  | \__/  \/  |___ |__/ 
                                                                                                                                                                              
# DB RUN (DOT IT FIRST !!!!!!!!!!)

# db_cur.execute("SELECT * FROM Run")
# db_trucs = db_cur.fetchall()
# for db_truc in db_trucs:
	# runID = db_truc['runID']
	# runPath = db_truc['runPath']
	# if not os.path.exists(runPath):
		# print " \t %s NOT FOUND" % runPath
		# path_found = glob.glob('%s/*/%s' % (os.path.dirname(runPath),runID))
		# if path_found:
			# print "is now -> %s \n" % path_found[0]
			# db_cur.execute("UPDATE Run SET runPath='%s' WHERE runID='%s'" % (path_found[0],runID))
		# else:
			# path_found = glob.glob('%s/*/*/%s' % (os.path.dirname(runPath),runID))
			# if path_found:
				# print "is now -> %s \n" % path_found[0]
				# db_cur.execute("UPDATE Run SET runPath='%s' WHERE runID='%s'" % (path_found[0],runID))
			# else:
				# print "      -> ****** STILL NOT FOUND ****** \n"
# db_con.commit()
# db_con.close()

# # DB ANALYSIS BAM PATH
# count_correct = 0
# count_not_found = 0

# runid2path = {}
# db_cur.execute("SELECT * FROM Run")
# db_trucs = db_cur.fetchall()
# for db_truc in db_trucs:
	# runid2path[db_truc['runID']] = db_truc['runPath']

# db_cur.execute("SELECT * FROM Analysis")
# db_trucs = db_cur.fetchall()
# for db_truc in db_trucs:
	# analysisID = db_truc['analysisID']
	# bamPath = db_truc['bamPath']
	# bamName = os.path.basename(bamPath)
	# sampleFolder = os.path.basename(os.path.dirname(bamPath))
	# try:
		# corrected_runPath = runid2path[db_truc['run']]
	# except:
		# print "\n***********run %s not found\n" % db_truc['run']
		# continue
	# if not os.path.exists(bamPath):
		# print " \t %s NOT FOUND" % bamPath
		# # path_found = glob.glob('%s/%s*/%s' % (runPath,analysisID))
		# bam = '%s/%s/%s' % (corrected_runPath,sampleFolder,bamName)
		# if os.path.exists(bam):
			# print "is now -> %s \n" % bam
			# db_cur.execute("UPDATE Analysis SET bamPath='%s' WHERE analysisID='%s'" % (bam,analysisID))
			# count_correct+=1
		# else:
			# print "      -> ****** STILL NOT FOUND ****** \n"
			# count_not_found+=1

# print "\n %s bamPath were corrected" % count_correct
# print "%s bamPath were still not found" % count_not_found
# db_con.commit()
# db_con.close()



### GET GENOMIC DESCRIPTION FOR FALSEPOSTIVES FILE
# z = open('/media/stuff/FalsePositives_updated.tsv','w')
# csv_writer = csv.writer(z,delimiter='\t')
# with open('/media/stuff/FalsePositives.tsv','r') as fp:
	# csv_reader = csv.reader(fp,delimiter='\t')
	# for line in csv_reader:
		# if line[0].startswith('#'):
			# continue
		# transcript_without_version = line[0]
		# transcript_description = line[1]
		# db_cur.execute("""SELECT DISTINCT chromosome,genomicDescription FROM Variant 
		# INNER JOIN VariantAnnotation ON VariantAnnotation.variant = Variant.variantID
		# WHERE transcript like '%s.%%' and transcriptDescription = '%s'""" % (transcript_without_version,transcript_description))
		# db_variant = db_cur.fetchone()
		# if not db_variant:
			# genomic_description = 'not found in db'
		# else:
			# genomic_description = '%s:%s' % (db_variant['chromosome'],db_variant['genomicDescription'])
		# csv_writer.writerow([transcript_without_version,transcript_description,genomic_description,line[3],line[4]])


### CORRECT GENOMIC DESCRIPTION WHEN chromosome is missing
db_cur.execute("SELECT variantID,chromosome,genomicDescription FROM Variant WHERE genomicDescription NOT LIKE 'chr%:%'")
db_variants = db_cur.fetchall()
for db_variant in db_variants:
	completeGenomicDescription = '%s:%s' % (db_variant['chromosome'],db_variant['genomicDescription'])
	db_cur.execute("UPDATE Variant SET genomicDescription='%s' WHERE variantID='%s'" % (completeGenomicDescription,db_variant['variantID']))
db_con.commit()
db_con.close()
