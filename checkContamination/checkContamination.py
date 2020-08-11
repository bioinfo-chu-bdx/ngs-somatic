#!/usr/bin/env python
from optparse import OptionParser
import subprocess
import openpyxl
import sqlite3
import pysam
import time
import json
import copy
import csv
import sys
import os

""" Script to analyse control contamination """


def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

def representsInt(s):
	try: 
		s = int(s)
		return s
	except ValueError:
		return s

def get_column2index(sheet): # renvoie dict col2index
	col2index = {}
	for i in range(1,sheet.max_column+1):
		col2index[sheet.cell(row=1,column=i).value] = i
	return col2index

def check_stderr(stderr_path,indent=0):
	with open(stderr_path,'r') as stderr:
		for line in stderr.readlines():
			line = line.replace('\n','')
			for i in range(indent):
				line = '\t%s' % line
			print line

############################################################################################
parser = OptionParser()
parser.add_option('-i', '--run-folder', help="Run folder", dest='run_folder')
(options, args) = parser.parse_args()

control_names = ['H2O','H20','NTC'] # liste des noms possibles pour les temoins negatifs

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

if os.path.isfile('%s/barcodes.json' % options.run_folder):
	with open(options.run_folder+'/barcodes.json', 'r') as g:
		barcodes_json = json.load(g)
else:
	print "error : barcodes.json not found in run folder"
	exit()

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

checkconta_folder = '%s/_checkContamination' % options.run_folder
if not os.path.isdir(checkconta_folder):
	subprocess.call(['mkdir', checkconta_folder])

# FIND CONTROL SAMPLES TO CHECK
control_barcode = {}
other_barcode = {}
for barcode in barcodes_json:
	sample = barcodes_json[barcode]['sample']
	sample_id = barcodes_json[barcode]['sample_id']
	sample_folder = '%s/%s' % (options.run_folder,sample)
	lib = barcodes_json[barcode]['barcode_description']
	target = barcodes_json[barcode]['target_region_filepath'].split('/')[-1]
	for rt in global_param['run_type']:
		if global_param['run_type'][rt]['target_bed'].split('/')[-1] == target:
			run_type = rt
			break
	checkConta_read_len = global_param['run_type'][run_type]['checkContamination_read_len']
	bam = '%s/%s/%s_%s.bam' % (options.run_folder,sample,sample,barcode)
	finalReport = '%s/%s/%s_%s.finalReport.xlsx' % (options.run_folder,sample,sample,barcode)

	if '-checkContamination' in sample:
		print "- note : sample %s is already a filtered one and will not be used" % sample
		continue
	if not os.path.isfile(bamfile):
		print "(warning : %s not found)" % bamfile
		continue

	iscontrol = False
	for control_name in control_names:
		if control_name in sample.upper():
			iscontrol = True
	if iscontrol:
		control_barcode[barcode] = {'sample':sample,'sample_id':sample_id,'sample_folder':sample_folder,'lib':lib,'target':target,'bam':bam,'finalreport':finalreport,'checkConta_read_len':checkConta_read_len}
	else:
		other_barcode[barcode] = {'sample':sample,'sample_id':sample_id,'sample_folder':sample_folder,'lib':lib,'target':target,'bam':bam,'finalreport':finalreport}

###  __   __   __   __   ___  __   __      ___       __           __   __       ___  __   __       
### |__) |__) /  \ /  ` |__  /__` /__`    |__   /\  /  ` |__|    /  ` /  \ |\ |  |  |__) /  \ |    
### |    |  \ \__/ \__, |___ .__/ .__/    |___ /~~\ \__, |  |    \__, \__/ | \|  |  |  \ \__/ |___ 

for barcode in control_barcode:

	###################
	### TARGET DATA ### 
	###################

	gene2transcript = {}
	db_cur.execute("SELECT * FROM Gene")
	db_genes = db_cur.fetchall()
	for gene in db_genes:
		if gene['geneID'] not in gene2transcript:
			gene2transcript[gene['geneID']] = gene['transcript']
	db_cur.execute("SELECT * FROM TargetedRegion WHERE panel = %s" % control_barcode[barcode]['target'])
	db_regions = db_cur.fetchall()
	amplicon_coverage = {}
	for region in db_regions:
		amplicon_coverage[region['targetedRegionName']] = {}
		amplicon_coverage[region['targetedRegionName']]['chr'] = region['chromosome']
		amplicon_coverage[region['targetedRegionName']]['start'] = region['start']
		amplicon_coverage[region['targetedRegionName']]['stop'] = region['stop']
		amplicon_coverage[region['targetedRegionName']]['gene_id'] = '%s_%s' % (region['gene'],region['details'])
		amplicon_coverage[region['targetedRegionName']]['transcript'] = gene2transcript[region['targetedRegionName']]['gene_id']]
		amplicon_coverage[region['targetedRegionName']]['contamination count'] = 0

	######################
	### BAM FILTERING  ###
	######################

	# CREATING FILTERED BAM FOLDER
	checkConta_bam_folder = '%s-checkContamination' % control_barcode[barcode]['sample_folder']
	checkConta_barcode = '%s-checkContamination' % barcode
	checkConta_sample = '%s-checkContamination' % control_barcode[barcode]['sample']
	checkConta_sample_id = '%s-checkContamination' % control_barcode[barcode]['sample_id']

	if not os.path.isdir(checkConta_bam_folder):
		print "- creating filtered bam folder..."
		subprocess.call(['mkdir',checkConta_bam_folder])

	# CREATING FILTERED BAM ENTRY IN BARCODES JSON
	if checkConta_barcode not in barcodes_json.keys():
		print "- adding filtered bam in barcodes_json..."
		barcodes_json[checkConta_barcode] = {}
		for key in barcodes_json[barcode].keys():
			barcodes_json[checkConta_barcode][key] = barcodes_json[barcode][key]
		barcodes_json[checkConta_barcode]['sample'] = checkConta_sample
		barcodes_json[checkConta_barcode]['sample_id'] = checkConta_sample_id
		json_text = json.dumps(barcodes_json, indent=4)
		bc_json = open(run_folder+'/barcodes.json','w')
		bc_json.write(json_text)
		bc_json.close()

	bamfile = pysam.Samfile(control_barcode[barcode]['bam'],'rb')
	checkConta_bam = '%s/%s_%s.bam' % (checkConta_bam_folder,checkConta_sample,checkConta_barcode)
	checkConta_bamfile = pysam.Samfile(checkConta_bam, 'wb', template=bamfile)

	### CREATE NEW BAM WITH READS > XX bases
	print "- bam filtering..."
	for read in bamfile.fetch():
		if (len(read.query) >= int(control_barcode[barcode]['checkConta_read_len'])): 	# aligned portion of the read, exclude soft-clipped bases
			checkConta_bamfile.write(read)
	bamfile.close()
	checkConta_bamfile.close()

	# SORT NEW BAM AND CREATE BAI
	print "- samtools sort..."
	subprocess.call(['samtools','sort',checkConta_bam,'-o',checkConta_bam])
	print "- samtools index..."
	subprocess.call(['samtools','index',checkConta_bam])

	########################
	### NEW BAM ANALYSIS ###
	########################

	print "- filtered bam complete analysis..."
	cmd = subprocess.call(['python','%s/run_analysis_illumina.py' % pipeline_folder,'--bam',checkConta_bam],stdout=open('%s/checkConta_bam_analysis.stdout.txt' % checkConta_bam_folder,'w'),stderr=open('%s/checkConta_bam_analysis.stderr.txt' % checkConta_bam_folder,'w'))
	check_stderr('%s/checkConta_bam_analysis.stderr.txt' % checkConta_bam_folder, indent=1)

	##################################################################
	### COMPARAISON AVEC COUVERTURE / VARIANTS DES AUTRES PATIENTS ###
	##################################################################

	print "- comparing with other samples..."
	# creation dic avec variants de l'annotation du check conta
	variants_in_bam_checkConta = {}
	bam_checkConta_finalreport = openpyxl.load_workbook('%s/%s_%s.finalReport.xlsx' % (checkConta_bam_folder,checkConta_sample,checkConta_barcode))
	bam_checkConta_covSheet = bam_checkConta_finalreport.get_sheet_by_name('Amplicon Coverage')
	bam_checkConta_annoSheet = bam_checkConta_finalreport.get_sheet_by_name('Annotation')
	bam_checkConta_covcol = get_column2index(bam_checkConta_covSheet)
	bam_checkConta_anncol = get_column2index(bam_checkConta_annoSheet)
	for i in range(2,bam_checkConta_covSheet.max_row+1):
		amplicon_coverage[bam_checkConta_covSheet.cell(row=i,column=bam_checkConta_covcol['region_id']).value]['contamination count'] = int(bam_checkConta_covSheet.cell(row=i,column=bam_checkConta_covcol['total_reads']).value)
	for i in range(2,bam_checkConta_annoSheet.max_row+1):
		if bam_checkConta_annoSheet.cell(row=i,column=bam_checkConta_anncol['Transcript']).value:
			transcript = bam_checkConta_annoSheet.cell(row=i,column=bam_checkConta_anncol['Transcript']).value.split('.')[0]
			c_nomen = bam_checkConta_annoSheet.cell(row=i,column=bam_checkConta_anncol['c.']).value
			variants_in_bam_checkConta[(transcript,c_nomen)] = []

	# recuperation de la liste des patients a comparer
	sample2compare = []
	for other_barcode in barcodes_json:
		other_sample = barcodes_json[other_barcode]['sample']
		other_bed = barcodes_json[other_barcode]['target_region_filepath'].split('/unmerged/detail/')[-1]
		other_lib = barcodes_json[other_barcode]['barcode_description']
		if ((other_lib == library) and (other_bed == target_bed)) and (other_barcode != barcode) and (other_barcode != checkConta_barcode):
			sample2compare.append((other_sample,other_barcode))
			finalReport_path = '%s/%s/%s_%s.finalReport.xlsx' % (run_folder,other_sample,other_sample,other_barcode)
			if os.path.isfile(finalReport_path):
				finalReport = openpyxl.load_workbook(finalReport_path)
				covSheet = finalReport.get_sheet_by_name('Amplicon Coverage')
				annoSheet = finalReport.get_sheet_by_name('Annotation')		
				covcol = get_column2index(covSheet)
				anncol = get_column2index(annoSheet)		
				for i in range(2,covSheet.max_row+1):
					amplicon = covSheet.cell(row=i,column=covcol['region_id']).value
					totalReads = covSheet.cell(row=i,column=covcol['total_reads']).value
					amplicon_coverage[amplicon][other_sample] = int(totalReads)
				for i in range(2,annoSheet.max_row+1):
					if annoSheet.cell(row=i,column=anncol['Transcript']).value:
						transcript = annoSheet.cell(row=i,column=anncol['Transcript']).value.split('.')[0]
						c_nomen = annoSheet.cell(row=i,column=anncol['c.']).value
						if (transcript,c_nomen) in variants_in_bam_checkConta.keys():
							variants_in_bam_checkConta[(transcript,c_nomen)].append('%s_%s' % (other_sample,other_barcode))

	###############################
	#### ecriture des resultats ###
	###############################

	print "- writing results..."
	conta_report = openpyxl.Workbook()
	ampliconSheet = conta_report.create_sheet(title='Amplicons Contamination')
	annotationSheet = conta_report.create_sheet(title='Annotation')
	try:
		del conta_report['Sheet']
	except:
		pass

	# AMPLICON SHEET
	header = ['Amplicon','Chr','Start','End','Gene_id', 'Transcript', '%s reads > %s b' % (sample,options.read_len)]
	for s in sample2compare:
		other_sample = s[0]
		header.append(other_sample)
	for i in range(len(header)):
		ampliconSheet.cell(row=1,column=i+1).value = header[i]
		ampliconSheet.cell(row=1,column=i+1).font = openpyxl.styles.Font(name='Calibri', size=11, bold=True)
		ampliconSheet.cell(row=1,column=i+1).border = openpyxl.styles.Border(left=openpyxl.styles.Side(style='thin'),right=openpyxl.styles.Side(style='thin'), top=openpyxl.styles.Side(style='thin'),bottom=openpyxl.styles.Side(style='thin'))

	to_write = []
	for amplicon in amplicon_coverage.keys():
		line = [amplicon,amplicon_coverage[amplicon]['chr'],amplicon_coverage[amplicon]['start'],amplicon_coverage[amplicon]['stop'],amplicon_coverage[amplicon]['gene_id'],amplicon_coverage[amplicon]['transcript'],amplicon_coverage[amplicon]['contamination count']]
		for s in sample2compare:
			try:
				line.append(amplicon_coverage[amplicon][s[0]])
			except:
				line.append('?')
		to_write.append(line)
	to_write.sort(key=lambda x: x[6]) # classement par ordre decroissant de nombre de reads > read_len
	to_write.reverse()

	for i in range(len(to_write)):
		for j in range(len(to_write[i])):
			ampliconSheet.cell(row=i+2,column=j+1).font = openpyxl.styles.Font(name='Calibri', size=11)
			if j <= 6:
				ampliconSheet.cell(row=i+2,column=j+1).value = representsInt(to_write[i][j])
				if j == 6 and int(to_write[i][6]) >= 100:
					ampliconSheet.cell(row=i+2,column=j+1).font = openpyxl.styles.Font(name='Calibri', size=11, color='ff0000')
			elif j > 6:
				try:
					if int(to_write[i][j]) < 5*(int(to_write[i][6])):
						ampliconSheet.cell(row=i+2,column=j+1).value = representsInt(to_write[i][j])
						ampliconSheet.cell(row=i+2,column=j+1).fill = openpyxl.styles.PatternFill(fill_type='solid',start_color='d28e8e')
					else:
						ampliconSheet.cell(row=i+2,column=j+1).value = representsInt(to_write[i][j])
				except:
					ampliconSheet.cell(row=i+2,column=j+1).fill = openpyxl.styles.PatternFill(fill_type='solid',start_color='d28e8e')

	# ANNOTATION SHEET
	annotationSheet.cell(row=1,column=1).value = 'Library Search'
	annotationSheet.cell(row=1,column=1).font = openpyxl.styles.Font(name='Calibri', size=11, bold=True)
	annotationSheet.cell(row=1,column=1).border = openpyxl.styles.Border(left=openpyxl.styles.Side(style='thin'),right=openpyxl.styles.Side(style='thin'), top=openpyxl.styles.Side(style='thin'),bottom=openpyxl.styles.Side(style='thin'))
	for j in range(1,bam_checkConta_annoSheet.max_column+1):
		annotationSheet.cell(row=1,column=j+1).value = bam_checkConta_annoSheet.cell(row=1,column=j).value
		annotationSheet.cell(row=1,column=j+1).font = openpyxl.styles.Font(name='Calibri', size=11, bold=True)
		annotationSheet.cell(row=1,column=j+1).border = openpyxl.styles.Border(left=openpyxl.styles.Side(style='thin'),right=openpyxl.styles.Side(style='thin'), top=openpyxl.styles.Side(style='thin'),bottom=openpyxl.styles.Side(style='thin'))

	for i in range(2,bam_checkConta_annoSheet.max_row+1):
		if bam_checkConta_annoSheet.cell(row=i,column=bam_checkConta_anncol['Transcript']).value:
			if not bam_checkConta_annoSheet.cell(row=i,column=bam_checkConta_anncol['Commentaire']).value:
				annotationSheet.cell(row=i,column=bam_checkConta_anncol['Commentaire']+1).value = '.'
			transcript = bam_checkConta_annoSheet.cell(row=i,column=bam_checkConta_anncol['Transcript']).value.split('.')[0]
			c_nomen = bam_checkConta_annoSheet.cell(row=i,column=bam_checkConta_anncol['c.']).value
			if variants_in_bam_checkConta[(transcript,c_nomen)]: # list not empty = variant found somewhere
				foundstring = ','.join(variants_in_bam_checkConta[(transcript,c_nomen)])
				annotationSheet.cell(row=i,column=1).value = 'Variant found in ' + foundstring
				annotationSheet.cell(row=i,column=1).font = openpyxl.styles.Font(name='Calibri', size=11, color='ff0000')
			else:
				annotationSheet.cell(row=i,column=1).value = 'not found'
				annotationSheet.cell(row=i,column=1).font = openpyxl.styles.Font(name='Calibri', size=11, color='00b050')
		for j in range(1,bam_checkConta_annoSheet.max_column+1):
			if bam_checkConta_annoSheet.cell(row=i,column=j).value:
				annotationSheet.cell(row=i,column=j+1).value = representsInt(bam_checkConta_annoSheet.cell(row=i,column=j).value)
				annotationSheet.cell(row=i,column=j+1).font = openpyxl.styles.Font(name='Calibri', size=11)
				annotationSheet.cell(row=i,column=j+1).fill = copy.copy(bam_checkConta_annoSheet.cell(row=i,column=j).fill)

	# coloration rouge cosmic
	for i in range(2,bam_checkConta_annoSheet.max_row+1):
		if annotationSheet.cell(row=i,column=bam_checkConta_anncol['COSMIC']+1).value and annotationSheet.cell(row=i,column=bam_checkConta_anncol['COSMIC']+1).value != '.':
			annotationSheet.cell(row=i,column=bam_checkConta_anncol['COSMIC']+1).font = openpyxl.styles.Font(name='Calibri', size=11, color='ff0000')

	## LAYOUT
	ws = conta_report['Amplicons Contamination']
	ws.column_dimensions[ws['A1'].column].width = 18
	ws.column_dimensions[ws['B1'].column].width = 8
	ws.column_dimensions[ws['C1'].column].width = 12
	ws.column_dimensions[ws['D1'].column].width = 12
	ws.column_dimensions[ws['E1'].column].width = 15
	ws.column_dimensions[ws['F1'].column].width = 15
	ws = conta_report['Annotation']
	maxsize = 15
	dims = {}
	for row in ws.rows:
		for cell in row:
			if cell.value:
				try:
					text_size = len(str(cell.value)) + 2
					dims[cell.column] = min(max(dims.get(cell.column, 0), text_size),maxsize)
				except:
					pass
	for col, value in dims.items():
		ws.column_dimensions[col].width = value
	ws.column_dimensions[ws['J1'].column].width = 8
					
	conta_report.save(checkconta_folder+'/checkContamination_%s.xlsx'%sample)
	subprocess.call(['mv',checkConta_bam_folder,checkconta_folder])
