#!/usr/bin/env python
import sys
import json
import os
import csv
import urllib2

# THIS SCRIPT NEED THE FOLLOWING FILES (decompressed with GUNZIP) : 
# - ncbiRefSeqCurated 	$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ncbiRefSeqCurated.txt.gz %s/reference_files
# - knownToRefSeq 		$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownToRefSeq.txt.gz %s/reference_files
# - knownCanonical 		$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownCanonical.txt.gz %s/reference_files

# USAGE : python annotate_bed.py /path/to/target.bed
# RESULTS : path/to/target.annotated.bed
# BED SHOULD HAVE LINUX LINE ENDING

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



pipeline_folder =  os.environ['NGS_PIPELINE_BX_DIR']
knownToRefSeq_path = '%s/reference_files/knownToRefSeq.txt' % pipeline_folder
knownCanonical_path = '%s/reference_files/knownCanonical.txt' % pipeline_folder
refGene_path = '%s/reference_files/ncbiRefSeqCurated.txt' % pipeline_folder
favorite_nm_path = '%s/reference_files/annotate_bed/favorite_NM_list.tsv' % pipeline_folder



bed_path = sys.argv[1]
bed_reader = open(bed_path,'r')
bed_writer = open(bed_path.split('.bed')[0]+'.anno.bed','w')

knownToRefSeq_file = open(knownToRefSeq_path,'r')
knownToRefSeq_reader = csv.reader(knownToRefSeq_file,delimiter='\t')
knownCanonical_file = open(knownCanonical_path,'r')
knownCanonical_reader = csv.reader(knownCanonical_file,delimiter='\t')
refGene_file = open(refGene_path,'r')
refGene_reader = csv.reader(refGene_file,delimiter='\t')
favorite_nm_file = open(favorite_nm_path,'r')
favorite_nm_reader = csv.reader(favorite_nm_file,delimiter='\t')

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

# PARSE favorite NMs
favorite_NM = {}
for fline in favorite_nm_reader:
	favorite_NM[fline[0]] = fline[1]

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
	if gene in favorite_NM.keys():
		print "\t- using favorite transcript"
		if transcript != favorite_NM[gene]:
			transcript = favorite_NM[gene]
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
	
	print "\t- GENE=%s;DETAILS=%s;TRANSCRIPT=%s;STRAND=%s" % (gene,details,transcript,strand)
	line2write.append("GENE=%s;DETAILS=%s;TRANSCRIPT=%s;STRAND=%s" % (gene,details,transcript,strand))
	l2w = '\t'.join(line2write)
	l2w = l2w+'\n'
	bed_writer.write(l2w)
	
bed_reader.close()
