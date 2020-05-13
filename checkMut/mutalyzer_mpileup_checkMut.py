#!/usr/bin/env python
import sys
import pysam
import os
from optparse import OptionParser
import subprocess
import json
import urllib2
import glob
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

############################################################################################

def convert_cpos_to_genomic(nm,cpos,version):
	request_sucess = False
	for v in range(int(version),0,-1):
		try: 
			url_mutalyzer_conversion = "https://mutalyzer.nl/json/numberConversion?build=hg19;variant=%s.%s:%s" % (nm,v,cpos.replace('+','%2B').replace('-','%2D'))
			f = urllib2.urlopen(url_mutalyzer_conversion)
			d = json.loads(f.read())
			chrom = d[0].split('NC_')[-1].split('.')[0]
			if (int(chrom) == 23):
				chrom = 'chrX'
			elif (int(chrom) == 24):
				chrom = 'chrY'
			else:
				chrom = 'chr' + str(int(chrom))
			start = int(d[0].split('g.')[-1].split('_')[0].split('A')[0].split('T')[0].split('G')[0].split('C')[0].split('del')[0].split('dup')[0])
			stop = start
			if '_' in d[0].split('g.')[-1]:
				stop = int(d[0].split('g.')[-1].split('_')[1].split('A')[0].split('T')[0].split('G')[0].split('C')[0].split('del')[0].split('ins')[0].split('dup')[0])
			request_sucess = True
		except:
			continue
		break
	if not request_sucess:
		print "* Error : in convert_cpos_to_genomic() : problem fetching mutalyzer data. (check if mutalyzer services are down, 502 bad gateway etc)"
	return (chrom,start,stop,v)
	
def get_offset_dup_del(chrom,start,stop):	# determiner si sequence est situe dans homopol et calculer decalage necessaire
	left_offset = 0
	right_offset = 0
	base = pysam.faidx(ref,'%s:%s-%s' % (chrom,start,stop)).split('\n')[1]
	for i in range(len(base),100,len(base)):
		left_seq = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-i,stop-i)).split('\n')[1]
		if left_seq == base:
			left_offset = i
		else:
			break
	for i in range(len(base),100,len(base)):
		right_seq = pysam.faidx(ref,'%s:%s-%s' % (chrom,start+i,stop+i)).split('\n')[1]
		if right_seq == base:
			right_offset = i
		else:
			break
	return (left_offset,right_offset)
	
def reverse(seq):
	seq = seq.replace('A','W').replace('T','X').replace('G','Y').replace('C','Z')
	seq = seq.replace('W','T').replace('X','A').replace('Y','C').replace('Z','G')
	seq = seq[::-1]
	return seq
	
def printing(text,mytext):
	print text
	txtoutput.write(text + '\n')
	mytext = mytext + text + '\n'
	return mytext
	
# ===============================================================================
parser = OptionParser()
parser.add_option('-b', '--bam',		help='Single bam ', 								dest='bam')
parser.add_option('-f', '--run-folder', help='Run folder ', 								dest='run_folder')
parser.add_option('-n', '--nm', 		help='NM transcipt (ex : NM_005228)', 				dest='nm')
parser.add_option('-c', '--cpos', 		help='c position in HGVS (ex : c.38G>A)', 			dest='cpos')
parser.add_option('-g', '--gene', 		help='use gene instead of NM',			 			dest='gene', default=False)
parser.add_option('-t', '--run-type', 	help='filter sample by run type. ex:LAM,TP53', 		dest='runtype', default=False)
parser.add_option('-l', '--bam-list',	help='List of bam', 								dest='bam_list')
parser.add_option('-s', '--sub-folder',	help='sub-folder of results', 						dest='sub_folder')
parser.add_option('-m', '--min-depth',	help='sample min depth at position',				dest='min_depth', default=20)
(options, args) = parser.parse_args()

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))
with open(global_param['NM_data'], 'r') as nmdata:
	NM_version = json.load(nmdata)
ref = global_param['default_reference']

if options.bam:
	bamlist = [options.bam]
elif options.bam_list:
	bamlist = []
	with open(options.bam_list,'r') as bl:
		for line in bl:
			bam = line.replace('\n','')
			bamlist.append(bam)
	options.run_folder = '%s/scripts/tests' % pipeline_folder
else:
	bamlist = glob.glob(options.run_folder+'/*/*.bam')
	bamlist = [item for item in bamlist if not 'processed' in item]
	
if options.runtype and not options.bam_list:
	options.runtype.replace(' ','')
	options.runtype = options.runtype.split(',')
	filtered_bamlist = []
	with open(options.run_folder+'/barcodes.json', 'r') as g:
		barcodes_json = json.load(g)
	for bamfile in bamlist:
		barcode = 'IonXpress_' + bamfile.split('IonXpress_')[-1].split('.bam')[0]
		target = barcodes_json[barcode]['target_region_filepath'].split('/')[-1]
		for _run_type in global_param['run_type']:
			if global_param['run_type'][_run_type]['target_bed'].split('/')[-1] == target:
				if _run_type in options.runtype:
					 filtered_bamlist.append(bamfile)
				break
	bamlist = filtered_bamlist
if not bamlist:
	print "- CheckMut : runtype not found, exit"
	exit()

if options.run_folder:
	checkMut_folder = options.run_folder + '/_checkMut'
	if not os.path.isdir(checkMut_folder):
		subprocess.call(['mkdir',checkMut_folder])
	if options.sub_folder:
		checkMut_folder = '%s/%s' % (checkMut_folder,options.sub_folder)
		if not os.path.isdir(checkMut_folder):
			subprocess.call(['mkdir',checkMut_folder])

mytext = ""
# ===============================================================================

cpos = options.cpos
if options.gene:
	if options.nm:
		gene = options.gene
		nm = options.nm
	else:
		gene = options.gene
		for item in NM_version:
			if NM_version[item]['gene'] == gene:
				nm = item
				break
else:
	nm = options.nm
	gene = NM_version[nm]['gene']
	
txtoutput = open('%s/%s_%s.txt' % (checkMut_folder,gene,cpos),'w')
	
version = NM_version[nm]['version']
genomics_pos = convert_cpos_to_genomic(nm,cpos,version)
chrom = genomics_pos[0]
start = genomics_pos[1]
stop = genomics_pos[2]
version_position_converter = genomics_pos[3]
strand = NM_version[nm]['strand']

# nucleotide a gauche et a droite dans la sequence du variant a prendre a compte pour comparer
vsize = 1
# nucleotide a gauche et a droite dans la sequence de ref a prendre a compte pour comparer
wsize = 0
# decalage potientiel en cas de duplication (ou deletion dans homopol)
right_offset = 0 
left_offset = 0

if '>' in cpos:
	variant_type = 'snp'
	reference = pysam.faidx(ref,'%s:%s-%s' % (chrom,start,stop)).split('\n')[1]
	variant = cpos.split('>')[-1]
	if strand == 'reverse':
		#reference = reverse(reference)
		variant = reverse(variant)
	mytext = printing("* snp = %s>%s" % (reference,variant),mytext)
elif 'delins' in cpos:
	variant_type = 'delins'
	deleted_size = (stop-start)+1
	inserted_size = len(cpos.split('delins')[-1])
	#wsize = max(3,max(deleted_size,inserted_size))
	wsize = max(deleted_size,inserted_size)
	while True:
		reference_large_windows = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-wsize,stop+wsize)).split('\n')[1]
		left = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-vsize,start-1)).split('\n')[1]
		right = pysam.faidx(ref,'%s:%s-%s' % (chrom,stop+1,stop+vsize)).split('\n')[1]
		if strand == 'reverse':
			rev = reverse(cpos.split('delins')[-1])
			variant = left + rev + right
		else:
			variant = left + cpos.split('delins')[-1] + right
		if variant in reference_large_windows or len(set(list(variant)))==1:
			vsize += 1
			if vsize == wsize:
				wsize += 1
		else:
			break
	reference = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-vsize,stop+vsize)).split('\n')[1]
	mytext = printing("* delins seq = %s -> %s" % (reference,variant),mytext)
elif 'del' in cpos:
	if '_' in cpos:
		variant_type = 'delx'
		wsize = (stop-start)+1
		while True:
			reference_large_windows = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-wsize,stop+wsize)).split('\n')[1]
			reference = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-vsize,stop+vsize)).split('\n')[1]
			left = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-vsize,start-1)).split('\n')[1]
			right = pysam.faidx(ref,'%s:%s-%s' % (chrom,stop+1,stop+vsize)).split('\n')[1]
			variant = left + right
			if variant in reference_large_windows or len(set(list(variant)))==1:
				vsize += 1
				if vsize == wsize:
					wsize += 1
			else:
				break
		mytext = printing("* del seq = %s -> %s" % (reference,left + '-'*(stop-start+1) + right),mytext)
	else:
		offsets = get_offset_dup_del(chrom,start,stop)
		left_offset = offsets[0]
		right_offset = offsets[1]
		deleted_base = pysam.faidx(ref,'%s:%s-%s' % (chrom,start,stop)).split('\n')[1]
		if right_offset > 0 or left_offset > 0 :
			variant_type = 'delx'
			reference = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-left_offset,stop+right_offset)).split('\n')[1]
			left = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-left_offset-1,start-left_offset-1)).split('\n')[1]
			right = pysam.faidx(ref,'%s:%s-%s' % (chrom,stop+right_offset+1,stop+right_offset+1)).split('\n')[1]
			variant = left + deleted_base*(len(reference)-1) + right
			wsize = len(variant) # nouvelle taille variant avec repetitions inclus
			#mytext = printing("(left offset = %s , right offset = %s)" % (left_offset,right_offset),mytext)
			mytext = printing("* del seq = %s -> %s" % (left+reference+right,variant),mytext)
		else:
			variant_type = 'del1'
			variant = 'deletion'
elif 'ins' in cpos:
	variant_type = 'insertion'
	reference = pysam.faidx(ref,'%s:%s-%s' % (chrom,start,stop)).split('\n')[1]
	if strand == 'reverse':
		inserted = reverse(cpos.split('ins')[-1])
	else:
		inserted = cpos.split('ins')[-1]
	variant = reference[0] + inserted + reference[1]
	wsize = len(variant)
	while True:
		reference_large_windows = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-wsize,stop+wsize)).split('\n')[1]
		reference = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-vsize,stop+vsize)).split('\n')[1]
		left = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-vsize,start)).split('\n')[1]
		right = pysam.faidx(ref,'%s:%s-%s' % (chrom,stop,stop+vsize)).split('\n')[1]
		variant = left + inserted + right
		if variant in reference_large_windows or len(set(list(variant)))==1:
			vsize += 1
			if vsize == wsize:
				wsize += 1
		else:
			break
	mytext = printing("* ins seq = %s -> %s" % (reference,variant),mytext)
elif 'dup' in cpos:
	variant_type = 'dup'
	offsets = get_offset_dup_del(chrom,start,stop)
	left_offset = offsets[0]
	right_offset = offsets[1]
	#mytext = printing("(left offset = %s , right offset = %s)" % (left_offset,right_offset),mytext)
	reference = pysam.faidx(ref,'%s:%s-%s' % (chrom,start-left_offset,stop+right_offset)).split('\n')[1]
	variant = reference+pysam.faidx(ref,'%s:%s-%s' % (chrom,start,stop)).split('\n')[1]
	wsize = len(variant)
	mytext = printing("* dup seq = %s -> %s" % (reference,variant),mytext)
	
if options.runtype:
	if 'ABL1' in options.runtype:
		chrom = 'ABL1.E5E6.NM_005157'
		start = int(cpos.split('c.')[-1].split('_')[0].split('A')[0].split('T')[0].split('G')[0].split('C')[0].split('del')[0].split('ins')[0].split('dup')[0]) + 192
		stop = start
		if '_' in cpos:
			stop = int(cpos.split('c.')[-1].split('_')[1].split('A')[0].split('T')[0].split('G')[0].split('C')[0].split('del')[0].split('ins')[0].split('dup')[0]) + 192

mytext = printing("* %s:%s.%s:%s" % (gene,nm,version_position_converter,cpos),mytext)
mytext = printing("* %s at %s:%s-%s" % (variant_type,chrom,start,stop),mytext)
mytext = printing("* %s strand" % strand,mytext)
mytext = printing("* min depth = %s\n" % options.min_depth,mytext)
mytext = printing("= SAMPLES ====================== FREQ(%) ===== COUNT(FWD/REV) ===== DEPTH",mytext)

#####################################
# BAM reading and counting variants #
#####################################

sample_names = []
sample_freqs = []
for bam in sorted(bamlist):
	sample = bam.split('/')[-1].split('_IonXpress')[0]
	bamfile = pysam.AlignmentFile(bam,'rb')
	no_read = True # Dans le cas de l'eau, souvent aucun read
	var_count = 0
	depth = 0
	negCov = 0
	posCov = 0
	if variant_type == 'snp' or variant_type == 'del1':
		base_count = {'A':0,'T':0,'G':0,'C':0,'deletion':0,'insertion':0}
		for pileupcolumn in bamfile.pileup(contig=chrom, start=start-1, stop=stop, max_depth = 100000, truncate=True, min_base_quality=0, stepper='nofilter'):
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_refskip:
					continue
				no_read = False
				## Note : on ne prend pas en compte les deletions dans la profondeur pour le calcul de la frequence (comme sur IGV/Alamut)
				## SAUF si le variant cible est une deletion
				if pileupread.is_del:
					if variant == 'deletion':
						depth += 1
						var_count += 1
						if pileupread.alignment.is_reverse: #negative strand hit
							negCov += 1
						else:
							posCov += 1
				else:
					depth += 1
					base = pileupread.alignment.query_sequence[pileupread.query_position]
					if variant == base:
						var_count += 1
						if pileupread.alignment.is_reverse: #negative strand hit
							negCov += 1
						else:
							posCov += 1

	if variant_type == 'delins' or variant_type == 'delx' or variant_type == 'dup' or variant_type == 'insertion':
		for pileupcolumn in bamfile.pileup(contig=chrom, start=start-1-wsize, stop=start-wsize, max_depth = 100000, truncate=True, min_base_quality=0, stepper='nofilter'):
			depth = pileupcolumn.n
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_refskip or pileupread.is_del:
					continue
				no_read = False
				extract_seq = pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+(wsize*3)]
				#print extract_seq
				if variant in extract_seq:
					var_count += 1
					if pileupread.alignment.is_reverse: #negative strand hit
						negCov += 1
					else:
						posCov += 1
	# SAMPLE RESULTS
	if no_read:
		mytext = printing("- %0-30s no reads at position" % (sample[0:25]),mytext)
		continue
	elif depth < int(options.min_depth):
		mytext = printing("- %0-30s not enough depth (%s)" % (sample[0:25],depth),mytext)
		continue
	else:
		freq = float(var_count)/float(depth)
		freq_in_percent = freq*100.0
	cnfwrv = "%s(%s/%s)" % (var_count,posCov,negCov)
	mytext = printing("- %0-30s %0-13.2f %0-20s %s" % (sample[0:25],freq_in_percent,cnfwrv,depth),mytext)
	sample_names.append(sample[0:25])
	sample_freqs.append(freq_in_percent)

mytext = printing("=========================================================================",mytext)
###### STATS WITH MEDIAN ABSOLUTE DEVIATION instead of mean standard deviation
median = numpy.median(sample_freqs)
mytext = printing("* median = %.2f" % median,mytext)
absolute_deviations = [numpy.absolute(sample_freq - median) for sample_freq in sample_freqs]
mad = 1.4826 * numpy.median(absolute_deviations)
mytext = printing("* mad (median absolute deviation) = %.2f" % mad,mytext)
noise_red_threshold = median + mad
mytext = printing("* high noise (red zone) : %.2f - %.2f" % (0,noise_red_threshold),mytext)
noise_end_threshold = median + 3*mad
mytext = printing("* uncertitude (orange zone) : %.2f - %.2f" % (noise_red_threshold,noise_end_threshold),mytext)
mytext = printing("* high-confidence variants (green zone) : >%.2f" % noise_end_threshold,mytext)
mytext = printing("=========================================================================",mytext)

### MATPLOTLIB GRAPH ###

if options.run_folder:
	fig = plt.figure()
	N = len(sample_names)
	ind = numpy.arange(N)
	width = 0.4
	ax = fig.add_subplot(111)
	rects = ax.bar(ind, sample_freqs, width, facecolor='#9999ff')
	plt.title('%s:%s:%s' % (gene,nm,cpos))
	ax.set_xlabel('Samples',fontsize=12)
	ax.set_ylabel('Variant Freq (%)',fontsize=12)
	plt.xticks(ind,sample_names,fontsize=4)
	max_y = max(sample_freqs) #ax.set_xlim([0,max(tab_freq_del)])
	ax.set_ylim([0,min(100,max(2,max_y*1.5))])	#augmentation un peu de l'ordonne, trop courte
	if max_y < 0.5:
		ax.spines['left'].set_color('red')
		ax.xaxis.label.set_color('red')
		ax.tick_params(axis='y', colors='red')
	ax.axhspan(0, noise_red_threshold, facecolor='red', alpha=0.1, label = 'noise')#\n(0 - %.2f)' % noise_red_threshold)
	ax.axhspan(noise_red_threshold, noise_end_threshold, facecolor='orange', alpha=0.1, label = 'incertain') #\n(%.2f - %.2f)' % (noise_red_threshold,noise_end_threshold))
	ax.axhspan(noise_end_threshold, max(2,max_y*1.5), facecolor='green', alpha=0.1, label = 'high confidence')#\n(> %.2f)' % noise_end_threshold)
	ax.legend(prop={'size':8})
	#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

	monofont = {'fontname':'FreeMono'}
	plt.text(0.95, 0.5, mytext, fontsize=8, va='center', transform=plt.gcf().transFigure, **monofont) # ha='left', va='left'
	
	fig.autofmt_xdate()
	fig.savefig('%s/%s_%s.png' % (checkMut_folder,gene,cpos),bbox_inches='tight',dpi=300)
