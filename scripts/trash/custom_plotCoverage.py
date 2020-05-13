#!/usr/bin/env python
import os
import subprocess
import json
import math
import matplotlib
matplotlib.use('Agg')
from pylab import *
import glob
import csv
from optparse import OptionParser

'''Script to generate coverage plots and coverage files'''

def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] for i in range(wanted_parts) ]

############################################################################################
parser = OptionParser()
parser.add_option('-i', '--run-folder', help='Run folder ', dest='run_folder')
(options, args) = parser.parse_args()

with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)

if os.path.isfile(options.run_folder+'/barcodes.json'):
	with open(options.run_folder+'/barcodes.json', 'r') as g:
		barcodes_json = json.load(g)
else:
	print "error : barcodes.json not found in run folder"

control_names = ['H2O','H20','NTC'] # liste des noms possibles pour les temoins negatifs
pool_file_path = global_param['pool_file_path']#'/media/n06lbth/sauvegardes_pgm/SBT/Panel_en_cours/IAD119108_231/IAD119108_231_384WellPlateDataSheet.csv'
plotcov_dir = options.run_folder + '/_2018_reanalysis/_plotCoverage'
my_colors = ['black','gray','silver','brown','red','salmon','orange','tan','peachpuff','khaki','yellow','olive','green','turquoise','cyan','dodgerblue','navy','mediumpurple','darkviolet','magenta','lightpink']
lstyle_list = ['solid','dotted','dashed']
############################################################################################

bamlist = glob.glob(options.run_folder+'/*/*.bam')
bamlist = [item for item in bamlist if not 'processed' in item]

barcode2target = {}
barcode2runtype = {}
barcode2sample = {}
barcode2covfile = {}
barcode2bamfile = {}
for bamfile in bamlist:
	sample = bamfile.split('/')[-1].split('_IonXpress')[0]
	barcode = 'IonXpress_' + bamfile.split('IonXpress_')[-1].split('.bam')[0]
	barcode2sample[barcode] = sample
	barcode2bamfile[barcode] = bamfile
	target = barcodes_json[barcode]['target_region_filepath'].split('/')[-1]
	for _run_type in global_param['run_type']:
		if global_param['run_type'][_run_type]['target_bed'].split('/')[-1] == target:
			barcode2runtype[barcode] = _run_type
			barcode2target[barcode] = global_param['run_type'][_run_type]['target_bed']
			break
	covfile = '/DATA/work/results/%s_%s.bam/coverage/%s_%s.amplicon.cov.xls' % (sample.replace(' ','@'),barcode,sample,barcode)
	if os.path.isfile(covfile):
		barcode2covfile[barcode] = covfile

################################
# PROCESS POUR CHAQUE RUN TYPE #
################################
if not os.path.isdir(plotcov_dir):
	subprocess.call(['mkdir',plotcov_dir])

for runtype in list(set(barcode2runtype.values())):
	print "Processing plotCoverage for %s samples" % runtype
	plotcov_runtype_dir = '%s/%s' % (plotcov_dir,runtype)
	if not os.path.isdir(plotcov_runtype_dir):
		subprocess.call(['mkdir',plotcov_runtype_dir])
		
	ampliconlist = [] # genomic ordered list of panel amplicons, for ploting
	barcode_sample = [] # liste de tuples (barcode, sample) pour l'etape de ploting
	barcode2amplcov = {} # key = barcode, value = amplicons
	ampl2gene = {} # gene name of amplicon, for plotting
	target_bed = global_param['run_type'][runtype]['target_bed']
	with open(target_bed, 'r') as tb:
		tb_reader = csv.reader(tb,delimiter='\t')
		tb_reader.next()
		for line in tb_reader:
			ampliconlist.append(line[3])
			ampl2gene[line[3]] = line[7].split('GENE_ID=')[-1].split(';')[0]
		
	for barcode in barcode2covfile.keys():
		sample = barcode2sample[barcode]
		if barcode2runtype[barcode] != runtype: # check if barcode is in runtype
			continue
		iscontrol = False # check if barcode is not a control
		for control_name in control_names:
			if control_name in sample.upper():
				iscontrol = True
		if iscontrol:
			continue
			
		print " - analysing %s coverage" % sample
		barcode_sample.append((barcode,sample))
		barcode2amplcov[barcode] = {}
		for amplicon in ampliconlist:
			barcode2amplcov[barcode][amplicon] = ''
			
		# get total read for coverageAnalysis file
		with open(barcode2covfile[barcode],'r') as cov_file:
			cov_file_reader = csv.reader(cov_file, delimiter='\t')
			cov_file_reader.next() # header
			for line in cov_file_reader:
				amplicon_id = line[3]
				total_reads = line[9]
				if int(total_reads) != 0:
					barcode2amplcov[barcode][amplicon_id] = math.log10(int(total_reads))
				else:
					barcode2amplcov[barcode][amplicon_id] = 0
		
	################################
	# create amplicons cov graphes #
	################################
	print "\n- Generating plots..."
	
	ampl2pool = {}
	try:
		pool_file = open(pool_file_path,'r')
		pool_reader = csv.reader(pool_file)
		for pline in pool_reader:
			try:
				amp = pline[3]
				if amp in ampl2gene.keys():
					if '_1_A' in pline[0]:
						ampl2pool[amp] = 'Pool 1'
					elif '_2_A' in pline[0]:
						ampl2pool[amp] = 'Pool 2'
			except:
				pass
		pool_file.close()
	except:
		print "warning : reading pool file failed"
	
	# get ordered list of barcode (order = alphabetical corresponding sample name)
	sorted_barcode_sample = sorted(barcode_sample, key=lambda tup: tup[1])
	sorted_barcodelist = [item[0] for item in sorted_barcode_sample]
	
	# create combinaisons of (samples,amplicons) to plot
	subset_max_amplicon = 80
	subset_max_barcode = 10
	amplicon_subset_num = int(math.ceil((float(len(ampliconlist))/float(subset_max_amplicon))))
	barcode_subset_num = int(math.ceil((float(len(sorted_barcodelist))/float(subset_max_barcode))))
	amplicon_subsets = split_list(ampliconlist, wanted_parts=amplicon_subset_num)
	barcode_subsets = split_list(sorted_barcodelist, wanted_parts=barcode_subset_num)
	
	amplicon_subsets.insert(0,ampliconlist) # full ampliconlist if the first set
	barcode_subsets.insert(0,sorted_barcodelist) # full barcodelist is the first set

	for i in range(len(amplicon_subsets)):
		for j in range(len(barcode_subsets)):
			lstyle_index = 0
			nb = 0
			
			# Axe X (Amplicons)
			x_names = []
			for amp_name in amplicon_subsets[i]:
				x_names.append(amp_name)
			x_len = len(x_names)
			x_data = range(len(x_names))
			x = array(x_data)

			# Axe Y (Samples)
			for barcode in barcode_subsets[j]:
				y_data = []
				for amplicon in amplicon_subsets[i]:
					y_data.append(barcode2amplcov[barcode][amplicon])
				y_legend = barcode2sample[barcode]
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
			# Adding gene name in amplicon legend
			for k in range(len(x_names)):
				xname = x_names[k]
				if x_names[k] in ampl2gene:
					xname = xname + ' - ' + ampl2gene[x_names[k]]
				if x_names[k] in ampl2pool:
					xname = ampl2pool[x_names[k]]  + ' - ' +  xname
				x_names[k] = xname
			xticks(x_data,x_names,rotation=90,fontsize=20)
			yticks(fontsize=25)
			ylim(ymin=0,ymax=5)
			xlabel('AMPLICONS',fontsize=30)
			ylabel('TOTAL READS (log10)',fontsize=30)	
			axhline(y=3,color='black',linestyle='--',label='1000X')
			axhline(y=2.7,color='b',linestyle='--',label='500X')
			axhline(y=2.48,color='orange',linestyle='--',label='300X')
			axhline(y=2,color='r',linestyle='--',label='100X')
			lgd = legend(loc='center left',bbox_to_anchor=(1,0.5),fontsize=25)

			# Output
			if i == 0:
				astring = 'all_amplicons'
			else:
				astring = 'amplicons_S%s' % str(i+1)
			if j == 0:
				bstring = 'all_samples'
			else:
				bstring = 'samples_S%s' % str(j+1)
			fig.autofmt_xdate() # inclinaison des legendes
			fig.savefig('%s/%s_%s.png' % (plotcov_runtype_dir,astring,bstring), bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=80)
			close()
