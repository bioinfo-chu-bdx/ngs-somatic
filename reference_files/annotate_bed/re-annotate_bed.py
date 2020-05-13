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

# USAGE : python annotate_bed.py /path/to/target.bed*
# bed should have dos line ending

panels = [
	'%s/reference_files/Target_ColonLung_v10_IAD172906_231.bed' % os.environ['NGS_PIPELINE_BX_DIR'], 
	'%s/reference_files/LAM2018_IAD143291_236_Designed_with_NM.bed', % os.environ['NGS_PIPELINE_BX_DIR']
	'%s/reference_files/TP53.20140108.designed_with_NM.bed', % os.environ['NGS_PIPELINE_BX_DIR']
	'%s/reference_files/Lymphome_B_IAD119887_231_Designed.with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'],
	'%s/reference_files/Lymphome_T_IAD120574_238_Designed.with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'],
	'%s/reference_files/IAD83112_241_Designed.with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'],
	#'%s/reference_files/ABL1_NM_005157_Designed.with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'] # a faire manuellement
	'%s/reference_files/FLT3_IAD161204_182_Designed.with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'],
	'%s/reference_files/old_ref_files/IAD78219_237_Designed_with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'],
	'%s/reference_files/old_ref_files/IAD37093_Designed_with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'],
	'%s/reference_files/old_ref_files/IAD62716_182_Designed_with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'],
	'%s/reference_files/old_ref_files/IAD165023_231_Designed.with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'],
	'%s/reference_files/old_ref_files/IAD154118_231_Designed.with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'],
	'%s/reference_files/old_ref_files/IAD119108_231_Designed.with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'],
	'%s/reference_files/old_ref_files/IAD108862_231_Designed.with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'],
	'%s/reference_files/old_ref_files/IAD94971_233_Designed.with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR'],
	'%s/reference_files/old_ref_files/IAD72953_231_Designed.with_NM.bed' % os.environ['NGS_PIPELINE_BX_DIR']
]

for bed_path in panels:
	print bed_path
	bed_reader = open(bed_path,'r')
	bed_name = bed_path.split('/')[-1]
	bed_writer = open('%s/reference_files/reannoted/%s' % (os.environ['NGS_PIPELINE_BX_DIR'],bed_name),'w')

	knownToRefSeq_path = '%s/reference_files/knownToRefSeq.txt' % os.environ['NGS_PIPELINE_BX_DIR']
	knownCanonical_path = '%s/reference_files/knownCanonical.txt' % os.environ['NGS_PIPELINE_BX_DIR']
	refGene_path = '%s/reference_files/ncbiRefSeqCurated.txt' % os.environ['NGS_PIPELINE_BX_DIR']

	knownToRefSeq_file = open(knownToRefSeq_path,'r')
	knownToRefSeq_reader = csv.reader(knownToRefSeq_file,delimiter='\t')
	knownCanonical_file = open(knownCanonical_path,'r')
	knownCanonical_reader = csv.reader(knownCanonical_file,delimiter='\t')
	refGene_file = open(refGene_path,'r')
	refGene_reader = csv.reader(refGene_file,delimiter='\t')

	amplnum = 1

	header = bed_reader.next() # header
	if "\"\"" in header:
		header = header.replace('""','@@')
		header = header.replace('"','')
		header = header.replace('@@','"')
	bed_writer.write(header)
	
	for bed_line in bed_reader:
		bl = bed_line.replace('\n','').split('\t')
		ampl_chr = bl[0]
		ampl_start = bl[1]
		ampl_stop = bl[2]
		old_gene = bl[7].split('GENE_ID=')[-1].split(';')[0].split('_')[0]
		nm = bl[7].split('SUBMITTED_REGION=')[-1].split(';')[0]
		try:
			amplicon = bl[3]
		except:
			amplicon = 'ampl' + str(amplnum)
			amplnum = amplnum + 1
		line2write = [ampl_chr,ampl_start,ampl_stop,amplicon,'0','+','.']
		gene = ''
		transcript = ''
		strand = ''
		exon = ''
		intron = ''
		test_pos = (int(bl[1])+int(bl[2]))/2 # position milieu ampl

		# FIND POSSIBLE CANONICAL USCS TRANSCRIPTS
		canonicals = []
		knownCanonical_file.seek(0) # begin after header
		knownCanonical_reader.next()
		for knownCanonical_line in knownCanonical_reader:
			canonical_start = int(knownCanonical_line[1])
			canonical_stop = int(knownCanonical_line[2])
			if (knownCanonical_line[0] == ampl_chr) and (canonical_start <= test_pos and canonical_stop >= test_pos) :
				canonicals.append(knownCanonical_line[4])
		if canonicals == []:
			print "\t-no canonical found, elarging windows of search"
			knownCanonical_file.seek(0) # begin after header
			knownCanonical_reader.next()
			for knownCanonical_line in knownCanonical_reader:
				canonical_start = int(knownCanonical_line[1])-1000
				canonical_stop = int(knownCanonical_line[2])+1000
				if (knownCanonical_line[0] == ampl_chr) and (canonical_start <= test_pos and canonical_stop >= test_pos) :
					canonicals.append(knownCanonical_line[4])
		if canonicals == []:
			print "\t-WARNING : no canonicals found"
			continue
			
		# FIND CORRESPONDING transcript
		knownToRefSeq_file.seek(0) # begin after header
		knownToRefSeq_reader.next()
		for knownToRefSeq_line in knownToRefSeq_reader:
			if knownToRefSeq_line[0] in canonicals and knownToRefSeq_line[1].startswith('NM_'):
				transcript = knownToRefSeq_line[1]
				refGene_file.seek(0) # begin after header
				refGene_reader.next()
				transcript = nm
				for refgene_line in refGene_reader:
					if refgene_line[1].split('.')[0] == transcript and ('_' not in refgene_line[2]) :
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
						# TARGETED EXON ? METHODE 1
						for index in range(len(exonStarts)):
							if (int(exonStarts[index]) < int(ampl_start) < int(exonEnds[index])) or (int(exonStarts[index]) < int(ampl_stop) < int(exonEnds[index])):
								if strand == 'forward':
									exon = index + 1
								elif strand == 'reverse':
									exon = len(exonStarts)-index
								break
						# TARGETED EXON ? METHODE 2
						if exon == '': #previous method didnt work, maybe amplicon totally overlapping exon (start + end in introns)
							for index in range(len(exonStarts)):
								if (int(ampl_start) < int(exonStarts[index]) < int(ampl_stop)) or (int(ampl_start) < int(exonEnds[index]) < int(ampl_stop)):
									if strand == 'forward':
										exon = index + 1
									if strand == 'reverse':
										exon = len(exonStarts) - index
									break
						# INTRON ?
						if exon == '':
							for index in range(len(exonStarts)-1):
								if (int(exonEnds[index]) < int(ampl_start) < int(exonStarts[index+1])) or (int(exonEnds[index]) < int(ampl_stop) < int(exonStarts[index+1])):
									if strand == 'forward':
										intron = index + 1
									elif strand == 'reverse':
										intron = (len(exonStarts)-index)-1
									break
						break
				break
		if old_gene != gene:
			print "\t\t WARNING : different gene : %s / %s" % (old_gene,gene)
		if transcript == '':
			print "\t\t WARNING : canonicals not found in refGene : %s" % canonicals
		if exon != '':
			details = 'ex%s' % exon
		elif intron != '':
			details = 'in%s' % intron
		else:
			details = '?'
		line2write.append("GENE=%s;DETAILS=%s;TRANSCRIPT=%s;STRAND=%s" % (gene,details,transcript,strand))
		l2w = '\t'.join(line2write)
		l2w = l2w+'\n'
		bed_writer.write(l2w)
