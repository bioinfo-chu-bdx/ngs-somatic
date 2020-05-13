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
parser.add_option('-i', '--bam', 			help="Input bam file ", 			dest='bam')
parser.add_option('-l', '--read-len', 		help="Min read len for filtering ", dest='read_len')
(options, args) = parser.parse_args()

sample = options.bam.split('/')[-1].split('_IonXpress')[0]
if '-filtered' in sample:
	print "error : this sample (%s) is already a filtered one" % sample
	exit()
barcode = 'IonXpress_' + options.bam.split('IonXpress_')[-1].split('.bam')[0]
sample_folder = os.path.dirname(options.bam)
run_folder = os.path.dirname(sample_folder)
if run_folder.endswith('/'):
	run_name = os.path.basename(os.path.dirname(run_folder))
else:
	run_name = os.path.basename(run_folder)

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))
	
db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()	
	
if not os.path.isfile(run_folder+'/barcodes.json'):
	print "error : barcodes.json not found in run folder"
	exit()
	
checkconta_folder = '%s/_checkContamination' % run_folder
if not os.path.isdir(checkconta_folder):
	subprocess.call(['mkdir', checkconta_folder])
	
with open(run_folder+'/barcodes.json', 'r') as g:
	barcodes_json = json.load(g)
	target_bed = barcodes_json[barcode]['target_region_filepath'].split('/unmerged/detail/')[-1]
	library = barcodes_json[barcode]['barcode_description']
	for rt in global_param['run_type']:
		if global_param['run_type'][rt]['target_bed'].split('/')[-1] == target_bed:
			run_type = rt
			break

target_bed_unmerged = global_param['run_type'][run_type]['target_bed']
sample_id = barcodes_json[barcode]['sample_id']

# TARGET BED DATA
amplicon_coverage = {}
bedfile = open(target_bed_unmerged,'r') 
target_reader = csv.reader(bedfile, delimiter = '\t')
target_reader.next()
for row in target_reader:
	amplicon_coverage[row[3]] = {}
	amplicon_coverage[row[3]]['chr'] = row[0]
	amplicon_coverage[row[3]]['start'] = row[1]
	amplicon_coverage[row[3]]['stop'] = row[2]
	g = row[7].split('GENE=')[-1].split(';')[0]
	d = row[7].split('DETAILS=')[-1].split(';')[0]
	amplicon_coverage[row[3]]['gene_id'] = '%s_%s' % (g,d)
	amplicon_coverage[row[3]]['transcript'] = row[7].split('TRANSCRIPT=')[-1].split(';')[0]
	amplicon_coverage[row[3]]['contamination count'] = 0
bedfile.close()

######################
### BAM FILTERING  ###
######################

# CREATING FILTERED BAM FOLDER
filtered_bam_folder = '%s-filtered' % sample_folder
filtered_bam_barcode = '%s-checkConta' % barcode
filtered_bam_sample = '%s-filtered' % barcodes_json[barcode]['sample']
filtered_bam_sample_id = '%s-filtered' % barcodes_json[barcode]['sample_id']
if not os.path.isdir(filtered_bam_folder):
	print "- creating filtered bam folder..."
	subprocess.call(['mkdir',filtered_bam_folder])
	
# CREATING FILTERED BAM ENTRY IN BARCODES JSON
if filtered_bam_barcode not in barcodes_json.keys():
	print "- adding filtered bam in barcodes_json..."
	barcodes_json[filtered_bam_barcode] = {}
	for key in barcodes_json[barcode].keys():
		barcodes_json[filtered_bam_barcode][key] = barcodes_json[barcode][key]
	barcodes_json[filtered_bam_barcode]['sample_id'] = filtered_bam_sample_id
	json_text = json.dumps(barcodes_json, indent=4)
	bc_json = open(run_folder+'/barcodes.json','w')
	bc_json.write(json_text)
	bc_json.close()

bamfile = pysam.Samfile(options.bam,'rb')
filtered_bam = '%s/%s_%s.bam' % (filtered_bam_folder,filtered_bam_sample,filtered_bam_barcode)
filtered_bamfile = pysam.Samfile(filtered_bam, 'wb', template=bamfile)

### CREATE NEW BAM WITH READS > 100 nc
print "- bam filtering..."
for read in bamfile.fetch():
	if (len(read.query) >= int(options.read_len)): 	# aligned portion of the read, exclude soft-clipped bases
		filtered_bamfile.write(read)
bamfile.close()
filtered_bamfile.close()

# SORT NEW BAM AND CREATE BAI
print "- samtools sort..."
subprocess.call(['samtools','sort',filtered_bam,'-o',filtered_bam])
print "- samtools index..."
subprocess.call(['samtools','index',filtered_bam])

#########################
### NEW BAM ANALYSIS  ###
#########################
	
print "- filtered bam complete analysis..."
cmd = subprocess.call(['python','%s/run_analysis.py' % pipeline_folder,'--bam',filtered_bam],stdout=open('%s/filtered_bam_analysis.stdout.txt' % filtered_bam_folder,'w'),stderr=open('%s/filtered_bam_analysis.stderr.txt' % filtered_bam_folder,'w'))
check_stderr('%s/filtered_bam_analysis.stderr.txt' % filtered_bam_folder, indent=1)

##################################################################
### COMPARAISON AVEC COUVERTURE / VARIANTS DES AUTRES PATIENTS ###
##################################################################

print "- comparing with other samples..."
# creation dic avec variants de l'annotation du check conta
variants_in_bam_filtered = {}
bam_filtered_finalreport = openpyxl.load_workbook('%s/%s_%s.finalReport.xlsx' % (filtered_bam_folder,filtered_bam_sample,filtered_bam_barcode))
bam_filtered_covSheet = bam_filtered_finalreport.get_sheet_by_name('Amplicon Coverage')
bam_filtered_annoSheet = bam_filtered_finalreport.get_sheet_by_name('Annotation')
bam_filtered_covcol = get_column2index(bam_filtered_covSheet)
bam_filtered_anncol = get_column2index(bam_filtered_annoSheet)
for i in range(2,bam_filtered_covSheet.max_row+1):
	amplicon_coverage[bam_filtered_covSheet.cell(row=i,column=bam_filtered_covcol['region_id']).value]['contamination count'] = int(bam_filtered_covSheet.cell(row=i,column=bam_filtered_covcol['total_reads']).value)
for i in range(2,bam_filtered_annoSheet.max_row+1):
	if bam_filtered_annoSheet.cell(row=i,column=bam_filtered_anncol['Transcript']).value:
		transcript = bam_filtered_annoSheet.cell(row=i,column=bam_filtered_anncol['Transcript']).value.split('.')[0]
		c_nomen = bam_filtered_annoSheet.cell(row=i,column=bam_filtered_anncol['c.']).value
		variants_in_bam_filtered[(transcript,c_nomen)] = []

# recuperation de la liste des patients a comparer
sample2compare = []
for other_barcode in barcodes_json:
	other_sample = barcodes_json[other_barcode]['sample']
	other_bed = barcodes_json[other_barcode]['target_region_filepath'].split('/unmerged/detail/')[-1]
	other_lib = barcodes_json[other_barcode]['barcode_description']
	if ((other_lib == library) and (other_bed == target_bed)) and (other_barcode != barcode) and (other_barcode != filtered_bam_barcode):
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
					if (transcript,c_nomen) in variants_in_bam_filtered.keys():
						variants_in_bam_filtered[(transcript,c_nomen)].append('%s_%s' % (other_sample,other_barcode))

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
for j in range(1,bam_filtered_annoSheet.max_column+1):
	annotationSheet.cell(row=1,column=j+1).value = bam_filtered_annoSheet.cell(row=1,column=j).value
	annotationSheet.cell(row=1,column=j+1).font = openpyxl.styles.Font(name='Calibri', size=11, bold=True)
	annotationSheet.cell(row=1,column=j+1).border = openpyxl.styles.Border(left=openpyxl.styles.Side(style='thin'),right=openpyxl.styles.Side(style='thin'), top=openpyxl.styles.Side(style='thin'),bottom=openpyxl.styles.Side(style='thin'))

for i in range(2,bam_filtered_annoSheet.max_row+1):
	if bam_filtered_annoSheet.cell(row=i,column=bam_filtered_anncol['Transcript']).value:
		if not bam_filtered_annoSheet.cell(row=i,column=bam_filtered_anncol['Commentaire']).value:
			annotationSheet.cell(row=i,column=bam_filtered_anncol['Commentaire']+1).value = '.'
		transcript = bam_filtered_annoSheet.cell(row=i,column=bam_filtered_anncol['Transcript']).value.split('.')[0]
		c_nomen = bam_filtered_annoSheet.cell(row=i,column=bam_filtered_anncol['c.']).value
		if variants_in_bam_filtered[(transcript,c_nomen)]: # list not empty = variant found somewhere
			foundstring = ','.join(variants_in_bam_filtered[(transcript,c_nomen)])
			annotationSheet.cell(row=i,column=1).value = 'Variant found in ' + foundstring
			annotationSheet.cell(row=i,column=1).font = openpyxl.styles.Font(name='Calibri', size=11, color='ff0000')
		else:
			annotationSheet.cell(row=i,column=1).value = 'not found'
			annotationSheet.cell(row=i,column=1).font = openpyxl.styles.Font(name='Calibri', size=11, color='00b050')
	for j in range(1,bam_filtered_annoSheet.max_column+1):
		if bam_filtered_annoSheet.cell(row=i,column=j).value:
			annotationSheet.cell(row=i,column=j+1).value = representsInt(bam_filtered_annoSheet.cell(row=i,column=j).value)
			annotationSheet.cell(row=i,column=j+1).font = openpyxl.styles.Font(name='Calibri', size=11)
			annotationSheet.cell(row=i,column=j+1).fill = copy.copy(bam_filtered_annoSheet.cell(row=i,column=j).fill)

# coloration rouge cosmic
for i in range(2,bam_filtered_annoSheet.max_row+1):
	if annotationSheet.cell(row=i,column=bam_filtered_anncol['COSMIC']+1).value and annotationSheet.cell(row=i,column=bam_filtered_anncol['COSMIC']+1).value != '.':
		annotationSheet.cell(row=i,column=bam_filtered_anncol['COSMIC']+1).font = openpyxl.styles.Font(name='Calibri', size=11, color='ff0000')

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
subprocess.call(['mv',filtered_bam_folder,checkconta_folder])
