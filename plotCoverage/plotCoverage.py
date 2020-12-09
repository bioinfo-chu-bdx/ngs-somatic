#!/usr/bin/env python
import os
import sys
import csv
import glob
import json
import math
import numpy
import sqlite3
import zipfile
import subprocess
import matplotlib
matplotlib.use('Agg')
from pylab import *
from optparse import OptionParser


'''Script to generate coverage plots and coverage files'''

def split_list(alist, wanted_parts=1):
	length = len(alist)
	return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] for i in range(wanted_parts) ]


def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

############################################################################################
parser = OptionParser()
parser.add_option('-i', '--run-folder', help='Run folder ', dest='run_folder')
(options, args) = parser.parse_args()

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

if os.path.isfile(options.run_folder+'/barcodes.json'):
	with open(options.run_folder+'/barcodes.json', 'r') as g:
		barcodes_json = json.load(g)
else:
	print "\t -error : barcodes.json not found in run folder"

control_names = ['H2O','H20','NTC'] # liste des noms possibles pour les temoins negatifs
pool_file_path = global_param['pool_file_path'] #'/media/n06lbth/sauvegardes_pgm/SBT/Panel_en_cours/IAD119108_231/IAD119108_231_384WellPlateDataSheet.csv'
plotcov_dir = options.run_folder + '/_plotCoverage'
my_colors = ['black','gray','silver','brown','red','salmon','orange','tan','peachpuff','khaki','yellow','olive','green','turquoise','cyan','dodgerblue','navy','mediumpurple','darkviolet','magenta','lightpink']
lstyle_list = ['solid','dotted','dashed']
############################################################################################

for barcode in barcodes_json:
	cov_file = False
	if os.path.isfile('%s/%s/intermediate_files/coverage/%s_%s.amplicon.cov.xls' % (options.run_folder,barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode)):
		cov_file = open('%s/%s/intermediate_files/coverage/%s_%s.amplicon.cov.xls' % (options.run_folder,barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode),'r')
	elif os.path.isfile('%s/%s/intermediate_files/coverage/%s_%s.target.cov.xls' % (options.run_folder,barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode)):
		cov_file = open('%s/%s/intermediate_files/coverage/%s_%s.target.cov.xls' % (options.run_folder,barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode),'r')
	elif os.path.isfile('%s/%s/intermediate_files.zip' % (options.run_folder,barcodes_json[barcode]['sample'])):
		archive = zipfile.ZipFile('%s/%s/intermediate_files.zip' % (options.run_folder,barcodes_json[barcode]['sample']), 'r')
		if 'coverage/%s_%s.amplicon.cov.xls' % (barcodes_json[barcode]['sample'],barcode) in archive.namelist():
			cov_file = archive.open('coverage/%s_%s.amplicon.cov.xls' % (barcodes_json[barcode]['sample'],barcode))
		elif 'coverage/%s_%s.target.cov.xls' % (barcodes_json[barcode]['sample'],barcode) in archive.namelist():
			cov_file = archive.open('coverage/%s_%s.target.cov.xls' % (barcodes_json[barcode]['sample'],barcode))
	barcodes_json[barcode]['coverage_file'] = cov_file

#############################
# PROCESS POUR CHAQUE PANEL #
#############################
if not os.path.isdir(plotcov_dir):
	subprocess.call(['mkdir',plotcov_dir])

panels = list(set([barcodes_json[barcode]['panel'] for barcode in barcodes_json]))
for panel in panels:
	print "\t -processing plotCoverage for %s samples" % panel
	plotcov_panel_dir = '%s/%s' % (plotcov_dir,panel)
	if not os.path.isdir(plotcov_panel_dir):
		subprocess.call(['mkdir',plotcov_panel_dir])

	regionlist = [] # genomic ordered list of panel regions, for ploting
	barcode_sample_tuples = [] # liste de tuples (barcode, sample) pour l'etape de ploting

	unsorted_regionlist = []
	db_cur.execute("SELECT TargetedRegion.chromosome,start,targetedRegionName,gene,details FROM TargetedRegion INNER JOIN Panel ON TargetedRegion.panel = Panel.panelID INNER JOIN Transcript ON Transcript.transcriptID = TargetedRegion.transcript WHERE panel='%s' ORDER BY start" % panel)
	db_target_regions = db_cur.fetchall()

	regname2regdesc = {}
	for db_target_region in db_target_regions:
		region = '%s_%s_%s' % (db_target_region['targetedRegionName'],db_target_region['gene'],db_target_region['details'])
		regname2regdesc[db_target_region['targetedRegionName']] = region
		unsorted_regionlist.append((region,int(db_target_region['chromosome'].replace('chr','').replace('X','23').replace('Y','24')),int(db_target_region['start'])))
	regionlist = sorted(unsorted_regionlist, key = lambda x: (x[1], x[2]))
	regionlist = [item[0] for item in regionlist]

	for barcode in barcodes_json:
		if barcodes_json[barcode]['panel'] != panel: # check if barcode is in panel
			continue
		sample = barcodes_json[barcode]['sample']
		iscontrol = False # check if barcode is not a control
		for control_name in control_names:
			if control_name in sample.upper():
				iscontrol = True
		if iscontrol:
			continue

		print "\t -analysing %s coverage" % sample
		barcode_sample_tuples.append((barcode,sample))

		barcodes_json[barcode]['coverage'] = {}
		for region in regionlist:
			barcodes_json[barcode]['coverage'][region] = 0

		cov_file_reader = csv.reader(barcodes_json[barcode]['coverage_file'], delimiter='\t')
		cov_file_reader.next() # header
		for line in cov_file_reader:
			# region_name = line[3]
			total_reads = int(numpy.round(float(line[9])))
			region = regname2regdesc[line[3]]
			if total_reads > 0:
				barcodes_json[barcode]['coverage'][region] = math.log10(total_reads)
			else:
				barcodes_json[barcode]['coverage'][region] = 0

	##############################
	# create regions cov graphes #
	##############################

	print "\t -generating plots..."
	reg2pool = {}
	try:
		pool_file = open(pool_file_path,'r')
		pool_reader = csv.reader(pool_file)
		for pline in pool_reader:
			try:
				reg = pline[3]
				if reg in region2geneex.keys():
					if '_1_A' in pline[0]:
						reg2pool[reg] = 'Pool 1'
					elif '_2_A' in pline[0]:
						reg2pool[reg] = 'Pool 2'
			except:
				pass
		pool_file.close()
	except:
		print "\t -warning : reading pool file failed"

	# get ordered list of barcode (order = alphabetical corresponding sample name)
	sorted_barcode_sample_tuples = sorted(barcode_sample_tuples, key=lambda tup: tup[1])
	sorted_barcodelist = [item[0] for item in sorted_barcode_sample_tuples]

	# create combinaisons of (samples,regions) to plot
	subset_max_region = 80
	subset_max_barcode = 10
	region_subset_num = int(math.ceil((float(len(regionlist))/float(subset_max_region))))
	barcode_subset_num = int(math.ceil((float(len(sorted_barcodelist))/float(subset_max_barcode))))
	region_subsets = split_list(regionlist, wanted_parts=region_subset_num)
	barcode_subsets = split_list(sorted_barcodelist, wanted_parts=barcode_subset_num)

	region_subsets.insert(0,regionlist) # full regionlist if the first set
	barcode_subsets.insert(0,sorted_barcodelist) # full barcodelist is the first set

	for i in range(len(region_subsets)):
		for j in range(len(barcode_subsets)):
			lstyle_index = 0
			nb = 0
			
			# Axe X (regions)
			x_names = []
			for reg_name in region_subsets[i]:
				x_names.append(reg_name)
			x_len = len(x_names)
			x_data = range(len(x_names))
			x = array(x_data)

			# Axe Y (Samples)
			for barcode in barcode_subsets[j]:
				y_data = []
				for region in region_subsets[i]:
					y_data.append(barcodes_json[barcode]['coverage'][region])
				y_legend = barcodes_json[barcode]['sample']
				fig = figure(1,figsize=(65,30))
				y = array(y_data)

				# Color and style
				if nb >= len(my_colors):
					nb = 0
					lstyle_index = lstyle_index + 1
					if lstyle_index >= len(lstyle_list):
						lstyle_index = 0
				my_color = my_colors[nb]
				
				# Sample Plot
				plot(x,y,label=y_legend,color=my_color,linestyle=lstyle_list[lstyle_index],linewidth=3.5)
				nb = nb + 1 
			# Adding gene name in region legend
			for k in range(len(x_names)):
				xname = x_names[k]
				# if x_names[k] in region2geneex:
					# xname = xname + ' - ' + region2geneex[x_names[k]]
				if x_names[k] in reg2pool:
					xname = reg2pool[x_names[k]]  + ' - ' +  xname
				x_names[k] = xname
			xticks(x_data,x_names,rotation=90,fontsize=20)
			yticks(fontsize=25)
			ylim(ymin=0,ymax=5)
			xlabel('REGIONS',fontsize=30)
			ylabel('TOTAL READS (log10)',fontsize=30)
			axhline(y=3,color='black',linestyle='--',label='1000X')
			axhline(y=2.7,color='b',linestyle='--',label='500X')
			axhline(y=2.48,color='orange',linestyle='--',label='300X')
			axhline(y=2,color='r',linestyle='--',label='100X')
			lgd = legend(loc='center left',bbox_to_anchor=(1,0.5),fontsize=25)

			# Output
			if i == 0:
				astring = 'all_regions'
			else:
				astring = 'regions_S%02d' % (i+1)
			if j == 0:
				bstring = 'all_samples'
			else:
				bstring = 'samples_S%02d' % (j+1)
			fig.autofmt_xdate() # inclinaison des legendes
			fig.savefig('%s/%s_%s.png' % (plotcov_panel_dir,astring,bstring), bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=80)
			close()
