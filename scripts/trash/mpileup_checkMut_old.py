#!/usr/bin/env python
import pysam
from optparse import OptionParser
import json
import urllib2
import glob
import numpy

############################################################################################

def convert_cpos_to_genomic(nm,cpos,NM_version):
	for version in range(int(NM_version[options.nm]),0,-1):
		try: 
			url_mutalyzer_conversion = "https://mutalyzer.nl/json/numberConversion?build=hg19;variant=%s.%s:%s" % (options.nm,version,options.cpos.replace('+','%2B').replace('-','%2D'))
			#print url_mutalyzer_conversion
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
		except:
			continue
		break
	return (chrom,start,stop,version)
	
def get_offset_dup_del(chrom,start,stop):	# determiner si sequence est situe dans homopol et calculer decalage necessaire
	left_offset = 0
	right_offset = 0
	base = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start,stop)).split('\n')[1]
	for i in range(len(base),100,len(base)):
		left_seq = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start-i,stop-i)).split('\n')[1]
		if left_seq == base:
			left_offset = i
		else:
			break
	for i in range(len(base),100,len(base)):
		right_seq = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start+i,stop+i)).split('\n')[1]
		if right_seq == base:
			right_offset = i
		else:
			break
	return (left_offset,right_offset)

# ===============================================================================
parser = OptionParser()
parser.add_option('-b', '--bam',		help='Single bam ', 								dest='bam')
parser.add_option('-r', '--run-folder', help='Run folder ', 								dest='run_folder')
parser.add_option('-n', '--nm', 		help='NM transcipt (ex : NM_005228)', 				dest='nm')
parser.add_option('-c', '--cpos', 		help='c position in HGVS (ex : c.38G>A)', 			dest='cpos')
(options, args) = parser.parse_args()

with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)
with open(global_param['mutalyzer_nm_list'], 'r') as mlist:
	NM_version = json.load(mlist)

genomics_pos = convert_cpos_to_genomic(options.nm,options.cpos,NM_version)
chrom = genomics_pos[0]
start = genomics_pos[1]
stop = genomics_pos[2]
version = genomics_pos[3]

# nucleotide a gauche et a droite dans la sequence du variant a prendre a compte pour comparer
vsize = 2
# nucleotide a gauche et a droite dans la sequence de ref a prendre a compte pour comparer
wsize = 0
# decalage potientiel en cas de duplication (ou deletion dans homopol)
right_offset = 0 
left_offset = 0

if '>' in options.cpos:
	variant_type = 'snp'
	reference = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start,stop)).split('\n')[1]
	variant = options.cpos.split('>')[-1]
	print "* snp = %s>%s" % (reference,variant)
elif 'delins' in options.cpos:
	variant_type = 'delins'
	reference = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start-vsize,stop+vsize)).split('\n')[1]
	left = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start-vsize,start-1)).split('\n')[1]
	right = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,stop+1,stop+vsize)).split('\n')[1]
	variant = left + options.cpos.split('delins')[-1] + right
	wsize = stop - start
	print "* delins seq = %s -> %s" % (reference,variant)
elif 'del' in options.cpos:
	if '_' in options.cpos:
		variant_type = 'delx'
		wsize = stop - start
		reference_large_windows = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start-wsize,stop+wsize)).split('\n')[1]
		while True:
			reference = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start-vsize,stop+vsize)).split('\n')[1]
			left = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start-vsize,start-1)).split('\n')[1]
			right = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,stop+1,stop+vsize)).split('\n')[1]
			variant = left + right
			if variant in reference_large_windows or len(set(list(variant)))==1:
				vsize += 1
			else:
				break
		print "* del seq = %s -> %s" % (reference,left + '-'*(stop-start+1) + right)
	else:
		offsets = get_offset_dup_del(chrom,start,stop)
		left_offset = offsets[0]
		right_offset = offsets[1]
		if right_offset > 0 or left_offset > 0 :
			variant_type = 'delx'
			reference = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start-left_offset,stop+right_offset)).split('\n')[1]
			left = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start-left_offset-1,start-left_offset-1)).split('\n')[1]
			right = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,stop+right_offset+1,stop+right_offset+1)).split('\n')[1]
			variant = left + deleted_base*(len(reference)-1) + right
			wsize = len(variant)
			print "(left offset = %s , right offset = %s)" % (left_offset,right_offset)
			print "* del seq = %s -> %s" % (left+reference+right,variant)
		else:
			variant_type = 'del1'
			variant = 'deletion'
elif 'ins' in options.cpos:
	variant_type = 'insertion'
	reference = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start,stop)).split('\n')[1]
	variant = reference[0]+options.cpos.split('ins')[-1]+reference[1]
	wsize = len(variant)
	reference_large_windows = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start-wsize,stop+wsize)).split('\n')[1]
	while True:
		reference = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start-vsize,stop+vsize)).split('\n')[1]
		left = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start-vsize,start)).split('\n')[1]
		right = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,stop,stop+vsize)).split('\n')[1]
		variant = left + options.cpos.split('ins')[-1] + right
		if variant in reference_large_windows or len(set(list(variant)))==1:
			vsize += 1
		else:
			break
	print "* ins seq = %s -> %s" % (reference,variant)
elif 'dup' in options.cpos:
	variant_type = 'dup'
	offsets = get_offset_dup_del(chrom,start,stop)
	left_offset = offsets[0]
	right_offset = offsets[1]
	print "(left offset = %s , right offset = %s)" % (left_offset,right_offset)
	reference = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start-left_offset,stop+right_offset)).split('\n')[1]
	variant = reference+pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % (chrom,start,stop)).split('\n')[1]
	wsize = len(variant)
	print "* dup seq = %s -> %s" % (reference,variant)

mutname = '%s_%s' % (options.nm,options.cpos)
print "* %s.%s:%s" % (options.nm,version,options.cpos)
print "* %s at %s:%s-%s\n" % (variant_type,chrom,start,stop)

if options.bam:
	bamlist = [options.bam]
else:
	bamlist = glob.glob(options.run_folder+'/*/*.bam')
	bamlist = [item for item in bamlist if not 'processed' in item]

sample_freqs = {}
for bam in bamlist:
	sample = bam.split('/')[-1].split('_IonXpress')[0]
	bamfile = pysam.AlignmentFile(bam,'rb')
	no_read = True # Dans le cas de l'eau, souvent aucun read
	freqs = []
	if variant_type == 'snp' or variant_type == 'del1':
		base_count = {'A':0,'T':0,'G':0,'C':0,'deletion':0,'insertion':0}
		for pileupcolumn in bamfile.pileup(contig=chrom, start=start-1, stop=stop, max_depth = 100000, truncate=True, min_base_quality=0, stepper='nofilter'):
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_refskip:
					continue
				no_read = False
				if pileupread.is_del:
					base_count['deletion'] += 1
				else:
					if pileupread.indel :
						base_count['insertion'] += 1
					base = pileupread.alignment.query_sequence[pileupread.query_position]
					base_count[base] += 1
			depth = base_count['A'] + base_count['C'] + base_count['G'] + base_count['T']
			if variant_type == 'deletion':
				depth = depth + base_count['deletion']
			#print "- Depth at position %s = %s :" % (start,depth)
			#print "- A: %s, C: %s, G: %s, T: %s, del: %s, ins: %s" % (base_count['A'],base_count['C'],base_count['G'],base_count['T'],base_count['deletion'],base_count['insertion'])
			freq = float(base_count[variant])/float(depth)
			freqs.append(freq)

	if variant_type == 'delins' or variant_type == 'delx' or variant_type == 'dup' or variant_type == 'insertion':
		for pileupcolumn in bamfile.pileup(contig=chrom, start=start-wsize, stop=start-wsize+1, max_depth = 100000, truncate=True, min_base_quality=0, stepper='nofilter'):
			var_count = 0
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_refskip or pileupread.is_del:
					continue
				no_read = False
				extract_seq = pileupread.alignment.query_sequence[pileupread.query_position+1:pileupread.query_position+1+(wsize*3)]
				#print "variant = %s \n extract_seq = %s \n" % (variant,extract_seq)
				if variant in extract_seq:
					var_count += 1
			#print "var_count = %s" % var_count
			#print "depth = %s" % pileupcolumn.n
			freq = float(var_count)/float(pileupcolumn.n)
			freqs.append(freq)
			
	#if variant_type == 'insertion': # variant_type == 'dup' or 
		#for pileupcolumn in bamfile.pileup(contig=chrom, start=start-left_offset-2, stop=start-left_offset-1, max_depth = 100000, truncate=True, min_base_quality=0, stepper='nofilter'):
			#var_count = 0
			#for pileupread in pileupcolumn.pileups:
				#if pileupread.is_refskip or pileupread.is_del:
					#continue
				#no_read = False
				#extract_seq = pileupread.alignment.query_sequence[pileupread.query_position+1:pileupread.query_position+1+len(variant)] # et si insA dans homopol A?
				##print "variant = %s \n extract_seq = %s \n" % (variant,extract_seq)
				#if extract_seq == variant:
					#var_count += 1
			##print "var_count = %s" % var_count
			##print "depth = %s" % pileupcolumn.n
			#freq = float(var_count)/float(pileupcolumn.n)
			#freqs.append(freq)

	# SAMPLE RESULTS
	if no_read:
		freq_in_percent = 0.0
	else:
		freq_in_percent = numpy.around(numpy.mean(freqs)*100.0,1)
	print "- {0:40} {1} %".format(sample,str(freq_in_percent))
	
	sample_freqs[sample] = freq_in_percent

print "======================================="
mean_freq = numpy.mean(sample_freqs.values())
print "* mean_freq = %s" % mean_freq
standard_deviation_freq = numpy.std(sample_freqs.values()) # ,ddof=1
print "* standard deviation = %s " % standard_deviation_freq
threshold_for_outliers = mean_freq + 2*standard_deviation_freq

freq_noise = [x for x in sample_freqs.values() if (x < threshold_for_outliers)]
mean_noise = numpy.mean(freq_noise)
print "* mean_noise = %s" % mean_noise
standard_deviation_noise = numpy.std(freq_noise) # ,ddof=1
print "* standard deviation after outliers removal = %s " % standard_deviation_noise
noise_end_threshold = mean_noise + standard_deviation_noise
print "* noise end threshold = %s" % noise_end_threshold
high_confidence_threshold = mean_noise + 3*standard_deviation_noise
print "* true variant high-confidence threshold = %s" % high_confidence_threshold
print "======================================="

# MANQUE DELINS
# TESTER REVERSE
# autres details; freq del / depth..
