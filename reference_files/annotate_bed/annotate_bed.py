#!/usr/bin/env python
import json
import os
import csv
import sqlite3
from optparse import OptionParser

# THIS SCRIPT NEED THE FOLLOWING FILES (decompressed with GUNZIP) : 
# - ncbiRefSeqCurated 	$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ncbiRefSeqCurated.txt.gz %s/reference_files
# - knownToRefSeq 		$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownToRefSeq.txt.gz %s/reference_files
# - knownCanonical 		$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownCanonical.txt.gz %s/reference_files

# USAGE : python annotate_bed.py --i /path/to/target.bed --o /path/to/target.anno.bed
# !!! BED SHOULD HAVE LINUX LINE ENDING !!! (use dos2unix script)

def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

def find_refgene_data(transcript):
	gene = ''
	strand = ''
	exon = ''
	intron = ''
	if transcript in refgene.keys():
		gene = refgene[transcript]['gene']
		strand = refgene[transcript]['strand']
		exonStarts = refgene[transcript]['exonStarts']
		exonEnds = refgene[transcript]['exonEnds']
		
		# TARGETED EXON ? METHODE 1
		for index in range(len(exonStarts)):
			if (int(exonStarts[index]) < int(start) < int(exonEnds[index])) or (int(exonStarts[index]) < int(stop) < int(exonEnds[index])):
				if strand == 'forward':
					exon = index + 1
				elif strand == 'reverse':
					exon = len(exonStarts)-index
				break
		# TARGETED EXON ? METHODE 2
		if exon == '': #previous method didnt work, maybe amplicon totally overlapping exon (start + end in introns)
			for index in range(len(exonStarts)):
				if (int(start) < int(exonStarts[index]) < int(stop)) or (int(start) < int(exonEnds[index]) < int(stop)):
					if strand == 'forward':
						exon = index + 1
					if strand == 'reverse':
						exon = len(exonStarts) - index
					break
		# INTRON ?
		if exon == '':
			for index in range(len(exonStarts)-1):
				if (int(exonEnds[index]) < int(start) < int(exonStarts[index+1])) or (int(exonEnds[index]) < int(stop) < int(exonStarts[index+1])):
					if strand == 'forward':
						intron = index + 1
					elif strand == 'reverse':
						intron = (len(exonStarts)-index)-1
					break
		if exon != '':
			details = 'ex%s' % exon
		elif intron != '':
			details = 'in%s' % intron
		else:
			details = 'NONE'
		return gene,strand,exonStarts,exonEnds,details

	else :
		gene = 'NONE'
		strand = 'NONE'
		exonStarts = 'NONE'
		exonEnds = 'NONE'
		details = 'NONE'
		return gene,strand,exonStarts,exonEnds,details

parser = OptionParser()
parser.add_option('-i', '--input-bed',	help="path to bed",dest='bed')
parser.add_option('-o', '--output-bed',	help="output bed",dest='output')
parser.add_option('-c', '--copy-panel',	help="optionnal : use same transcripts as an existing panel in the DB",dest='copy_panel')
(options, args) = parser.parse_args()

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

favorite_transcript = {}
if options.copy_panel:
	db_cur.execute("""SELECT DISTINCT gene,transcriptID FROM Transcript 
	INNER JOIN TargetedRegion ON TargetedRegion.transcript = Transcript.transcriptID
	WHERE panel = '%s'
	ORDER BY Gene""" % options.copy_panel)
	db_transcripts = db_cur.fetchall()
	for db_transcript in db_transcripts:
		favorite_transcript[db_transcript['gene']] = db_transcript['transcriptID'].split('.')[0]

# favorite_transcript = {'ANKRD26':'NM_014915','ASXL1':'NM_015338','ASXL2':'NM_018263','BCOR':'NM_001123385','BCORL1':'NM_021946','CALR':'NM_004343',
# 'CBL':'NM_005188','CCND2':'NM_001759','CEBPA':'NM_004364','CSF3R':'NM_156039','CUX1':'NM_181552','DDX41':'NM_016222','DHX15':'NM_001358',
# 'DNMT3A':'NM_022552','ETNK1':'NM_018638','ETV6':'NM_001987','EZH2':'NM_004456','FLT3':'NM_004119','GATA1':'NM_002049','GATA2':'NM_032638',
# 'GNAS':'NM_000516','GNB1':'NM_002074','IDH1':'NM_005896','IDH2':'NM_002168','IKZF1':'NM_006060','JAK2':'NM_004972','KDM6A':'NM_021140',
# 'KIT':'NM_000222','KMT2C':'NM_170606','KRAS':'NM_033360','MPL':'NM_005373','MYC':'NM_001354870','NFE2':'NM_006163','NPM1':'NM_002520',
# 'NRAS':'NM_002524','PHF6':'NM_001015877','PPM1D':'NM_003620','PTEN':'NM_000314','PTPN11':'NM_002834','RAD21':'NM_006265','RIT1':'NM_006912',
# 'RUNX1':'NM_001754','SAMD9':'NM_017654','SAMD9L':'NM_001303497','SETBP1':'NM_015559','SF3B1':'NM_012433','SH2B3':'NM_005475','SMC1A':'NM_006306',
# 'SMC3':'NM_005445','SRSF2':'NM_003016','SRY':'NM_003140','STAG2':'NM_001042749','TERC':'NR_001566','TERT':'NM_198253','TET2':'NM_001127208',
# 'TP53':'NM_001126112','U2AF1':'NM_001025203','WT1':'NM_001198551','ZRSR2':'NM_005089'}

knownToRefSeq_path = global_param['knownToRefSeq']
knownCanonical_path = global_param['knownCanonical']
refGene_path = global_param['RefSeq']

bed_path = options.bed
bed_reader = open(bed_path,'r')
bed_writer = open(options.output,'w')

knownToRefSeq_file = open(knownToRefSeq_path,'r')
knownToRefSeq_reader = csv.reader(knownToRefSeq_file,delimiter='\t')
knownCanonical_file = open(knownCanonical_path,'r')
knownCanonical_reader = csv.reader(knownCanonical_file,delimiter='\t')
refGene_file = open(refGene_path,'r')
refGene_reader = csv.reader(refGene_file,delimiter='\t')

# PARSE CANONICAL
knownCanonical = []
for knownCanonical_line in knownCanonical_reader:
	canonical_chr = knownCanonical_line[0]
	canonical_start = int(knownCanonical_line[1])
	canonical_stop = int(knownCanonical_line[2])
	canonical = knownCanonical_line[4]
	knownCanonical.append([canonical_chr,canonical_start,canonical_stop,canonical])

# PARSE knownToRefSeq
knownToRefSeq = {}
for knownToRefSeq_line in knownToRefSeq_reader:
	if knownToRefSeq_line[1].startswith('NM_'):
		canonical = knownToRefSeq_line[0]
		transcript = knownToRefSeq_line[1]
		knownToRefSeq[canonical] = transcript

# PARSE refGene
refgene = {}
for refgene_line in refGene_reader:
	if refgene_line[1].startswith('NM_'):
		transcript = refgene_line[1].split('.')[0]
		gene = refgene_line[12]
		gene_strand = refgene_line[3]
		if gene_strand == '+':
			strand = 'forward'
		elif gene_strand == '-':
			strand = 'reverse'
		exonStarts = refgene_line[9].split(',')
		if '' in exonStarts:
			exonStarts.remove('')
		exonEnds = refgene_line[10].split(',')
		if '' in exonEnds:
			exonEnds.remove('')
		refgene[transcript] = {'gene':gene,'strand':strand,'exonStarts':exonStarts,'exonEnds':exonEnds}


# RE-WRITE THE HEADER IF EXISTANT
firstline = bed_reader.next()
if not firstline.startswith('chr'):
	header = firstline
	if "\"\"" in header:
		header = header.replace('""','@@')
		header = header.replace('"','')
		header = header.replace('@@','"')
	bed_writer.write(header)
else:
	bed_reader.seek(0)

regnum = 1

# ANNOTATE
for bed_line in bed_reader:
	bl = bed_line.replace('\n','').split('\t')
	chrom = bl[0]
	start = bl[1]
	stop = bl[2]
	print "- %s:%s-%s :" % (chrom,start,stop)
	try:
		region_id = bl[3]
	except:
		region_id = 'Region' + str(regnum)
		regnum = regnum + 1
	line2write = [chrom,start,stop,region_id,'0','+','.']
	gene = ''
	transcript = ''
	strand = ''
	exon = ''
	intron = ''
	test_pos = (int(bl[1])+int(bl[2]))/2 # position milieu ampl

	# FIND POSSIBLE CANONICAL USCS TRANSCRIPTS
	canonicals = []
	for kc in knownCanonical:
		if (kc[0] == chrom) and (kc[1] <= test_pos and kc[2] >= test_pos) :
			canonicals.append(kc[3])
	if canonicals == []:
		print "\t- elarging search windows..."
		z = 1
		while canonicals == [] and (z < 10):
			for kc in knownCanonical:
				if (kc[0] == chrom) and (kc[1]-(z*1000) <= test_pos and kc[2]+(z*1000) >= test_pos) :
					canonicals.append(kc[3])
			z +=1
	if canonicals == []:
		print "\t- WARNING : no canonicals found!"
		continue
		
	# FIND CORRESPONDING transcript
	for canonical in canonicals:
		if canonical in knownToRefSeq.keys():
			canonical_used = canonical
			transcript = knownToRefSeq[canonical_used]
			break
	
	# FIND REFGENE DATA
	gene, strand, exonStarts, exonEnds, details = find_refgene_data(transcript)
	if gene in favorite_transcript.keys():
		print "\t- using favorite transcript"
		if transcript != favorite_transcript[gene]:
			transcript = favorite_transcript[gene]
			gene, strand, exonStarts, exonEnds, details = find_refgene_data(transcript)
	else:
		if not canonicals:
			print "\t WARNING : canonicals not found between %s-%s" % (start,stop)
		else:
			if len(canonicals) > 1:
				print "\t WARNING : several canonicals found between %s-%s  : %s" % (start,stop,canonicals)
				print "\t Canonical used was %s" % canonical_used
			if transcript == '':
				print "\t WARNING : transcript not found in refGene for canonicals : %s" % canonicals

	line = "GENE=%s;DETAILS=%s;TRANSCRIPT=%s;STRAND=%s" % (gene,details,transcript,strand)
	if "Pool" in bl[5]:
		pool = bl[5].split("Pool=")[-1].split(";")[0]
		line = "%s;Pool=%s" % (line,pool)
	print "\t- %s" % line
	line2write.append(line)
	l2w = '\t'.join(line2write)
	l2w = l2w+'\n'
	bed_writer.write(l2w)
	
bed_reader.close()
