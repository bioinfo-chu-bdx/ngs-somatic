#!/usr/bin/env python
import os
import json
import glob
import pysam
import numpy
import sqlite3
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.variantmapper
import hgvs.exceptions
import hgvs.normalizer
import hgvs.exceptions
import hgvs.validator
import hgvs.parser


############################################################################################

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

def get_offset_dup_del(chrom,start,stop): # determiner si sequence est situe dans homopol et calculer decalage necessaire
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
parser.add_option('-f', '--run-folder', help='Run folder ',									dest='run_folder')
#parser.add_option('-n', '--transcript', help='NM transcipt (ex : NM_005228)',				dest='transcript')
parser.add_option('-c', '--cpos', 		help='c position in HGVS (ex : c.38G>A)',			dest='cpos')
parser.add_option('-g', '--gene', 		help='use gene instead of NM',			 			dest='gene', default=False)
parser.add_option('-t', '--run-type', 	help='filter sample by run type. ex:LAM,TP53',		dest='runtype', default=False)
parser.add_option('-s', '--sub-folder',	help='sub-folder of results',						dest='sub_folder')
parser.add_option('-m', '--min-depth',	help='sample min depth at position',				dest='min_depth', default=20)
parser.add_option('-z', '--min-sample', help='min samples to compare (use others runs)',	dest='min_sample', default=False)
parser.add_option('-p', '--use-processed', help='use processed bam',						dest='processed', default=False, action='store_true')
(options, args) = parser.parse_args()

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

ref = global_param['default_reference']

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

if not options.run_folder:
	print "- Error : need --run-folder arg"
	exit()

with open(options.run_folder+'/barcodes.json', 'r') as g:
	barcodes_json = json.load(g)

if options.runtype:
	options.runtype.replace(' ','')
	options.runtype = options.runtype.split(',')

bamlist = []
for barcode in barcodes_json:
	if options.processed:
		bamfile = '%s/%s/%s_%s.processed.bam' % (options.run_folder,barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode)
	else:
		bamfile = '%s/%s/%s_%s.bam' % (options.run_folder,barcodes_json[barcode]['sample'],barcodes_json[barcode]['sample'],barcode)
	if options.runtype:
		db_cur.execute("SELECT panelProject FROM Panel WHERE panelID='%s'" % barcodes_json[barcode]['panel'])
		db_panel = db_cur.fetchone()
		project = db_panel['panelProject']
		if project in options.runtype:
			if os.path.isfile(bamfile):
				bamlist.append(bamfile)
			else:
				print "(warning : %s not found)" % bamfile
	else:
		if os.path.isfile(bamfile):
			bamlist.append(bamfile)
		else:
			print "(warning : %s not found)" % bamfile

if not bamlist:
	print "- CheckMut : no BAM found to process, exit"
	exit()

added_bam = []
if options.min_sample:
	try:
		min_sample = int(options.min_sample)
		db_cur.execute("SELECT * FROM Analysis INNER JOIN Panel ON Panel.panelID = Analysis.panel INNER JOIN Sample ON Sample.sampleID = Analysis.sample WHERE panelProject IN %s AND isControl=0 ORDER BY analysisDate DESC" % str(options.runtype).replace('[','(').replace(']',')'))
		db_analyzes = db_cur.fetchall()
		for db_analysis in db_analyzes:
			bamPath = db_analysis['bamPath']
			if options.processed:
				bamPath = bamPath.replace('.bam','.processed.bam')
			if os.path.exists(bamPath) and bamPath not in bamlist:
				bamlist.append(bamPath)
				added_bam.append(bamPath)
			if len(bamlist) >= min_sample:
				break
	except Exception as e:
		print str(e)
		print "* Error in adding sample"

checkMut_folder = options.run_folder + '/_checkMut'
if not os.path.isdir(checkMut_folder):
	subprocess.call(['mkdir',checkMut_folder])
if options.sub_folder:
	checkMut_folder = '%s/%s' % (checkMut_folder,options.sub_folder)
	if not os.path.isdir(checkMut_folder):
		subprocess.call(['mkdir',checkMut_folder])

mytext = ""

hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
hn = hgvs.normalizer.Normalizer(hdp)
hv = hgvs.validator.Validator(hdp)
am37 = easyvariantmapper = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37')
# ===============================================================================

cpos = options.cpos
gene = options.gene

if options.run_folder.endswith('/'):
	run_name = options.run_folder.split('/')[-2]
else:
	run_name = options.run_folder.split('/')[-1]
db_cur.execute("""SELECT DISTINCT transcriptID, strand FROM Transcript
INNER JOIN TargetedRegion ON TargetedRegion.transcript = Transcript.transcriptID
INNER JOIN Panel ON Panel.panelID = TargetedRegion.panel
INNER JOIN Analysis On Analysis.panel = Panel.panelID
WHERE run = '%s' AND gene = '%s'""" % (run_name,options.gene))
db_transcript = db_cur.fetchone()
transcript = db_transcript['transcriptID']
strand = db_transcript['strand']

if os.path.isfile('%s/_checkMut.txt' % checkMut_folder):
	txtoutput = open('%s/_checkMut.txt' % checkMut_folder,'a')
else:
	txtoutput = open('%s/_checkMut.txt' % checkMut_folder,'w')
#txtoutput = open('%s/%s_%s.txt' % (checkMut_folder,gene,cpos),'w')

# c = hp.parse_hgvs_variant('%s.%s:%s' % (transcript, transcript_version, cpos))
c = hp.parse_hgvs_variant('%s:%s' % (transcript, cpos))
g = am37.c_to_g(c)
chrom = g.ac.split('.')[0][-2:] # ex : 04
chrom = int(chrom)
if chrom == 23:
	chrom = 'X'
elif chrom == 24:
	chrom = 'Y'
chrom = 'chr%s' % chrom
start = int(g.posedit.pos.start.base)
stop = int(g.posedit.pos.end.base)

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

# mytext = printing("* %s:%s.%s:%s" % (gene,transcript,transcript_version,cpos),mytext)
mytext = printing("* %s:%s:%s" % (gene,transcript,cpos),mytext)
mytext = printing("* %s at %s:%s-%s" % (variant_type,chrom,start,stop),mytext)
mytext = printing("* %s strand" % strand,mytext)
mytext = printing("* min depth = %s\n" % options.min_depth,mytext)
if added_bam:
	mytext = printing("* (+) means sample was added from other run",mytext)
mytext = printing("= SAMPLES ====================== FREQ(%) ===== COUNT(FWD/REV) ===== DEPTH",mytext)

#####################################
# BAM reading and counting variants #
#####################################

sample_names = []
sample_freqs = []
for bam in sorted(bamlist):
	sample = bam.split('/')[-1].split('_S')[0]
	if bam in added_bam:
		sample = '(+)%s' % sample
	try:
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
	except Exception as e:
		print "(Error processing sample %s : %s)" % (sample,e)

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
	plt.title('%s:%s:%s' % (gene,transcript,cpos))
	ax.set_xlabel('Samples',fontsize=12)
	ax.set_ylabel('Variant Freq (%)',fontsize=12)
	if N<=50:
		plt.xticks(ind,sample_names,fontsize=4)
	else:
		plt.xticks(ind,sample_names,fontsize=2)
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
	if N<=50:
		plt.text(0.95, 0.5, mytext, fontsize=8, va='center', transform=plt.gcf().transFigure, **monofont) # ha='left', va='left'
		fig.autofmt_xdate()
		fig.savefig('%s/%s_%s.png' % (checkMut_folder,gene,cpos.replace('>','_')),bbox_inches='tight',dpi=300)
	else:
		plt.text(0.95, 0.5, mytext, fontsize=4, va='center', transform=plt.gcf().transFigure, **monofont) # ha='left', va='left'
		fig.autofmt_xdate()
		fig.savefig('%s/%s_%s.png' % (checkMut_folder,gene,cpos.replace('>','_')),bbox_inches='tight',dpi=600)
