#!/usr/bin/env python
import sys
import json
import os
import csv
import urllib2

# TODO https://genome.ucsc.edu/cgi-bin/hgTables
# obtention fichier refseq : hg19 - Genes and Gene Predicitons - RefSeq Genes - refGene, puis supprimer les doublons NM avec excel
# obtention fichier knownToRefSeq : hg19 - Genes and Gene Predicitons - RefSeq Genes - knownToRefSeq
# obtention fichier knownCanonical : hg19 - Genes and Gene Predicitons - UCSC Genes - knownCanonical

bed_path = sys.argv[1]

knownToRefSeq_path = '/DATA/work/reference_files/Annotate_bed_script/knownToRefSeq.txt'
knownCanonical_path = '/DATA/work/reference_files/Annotate_bed_script/knownCanonical.txt'
refGene_path = '/DATA/work/reference_files/Annotate_bed_script/refGene_sans_doublons_NM.txt'

knownToRefSeq_file = open(knownToRefSeq_path,'r')
knownCanonical_file = open(knownCanonical_path,'r')
refGene_file = open(refGene_path,'r')
bed_file = open(bed_path,'r')

knownToRefSeq_reader = csv.reader(knownToRefSeq_file,delimiter='\t')
knownCanonical_reader = csv.reader(knownCanonical_file,delimiter='\t')
refGene_reader = csv.reader(refGene_file,delimiter='\t')
bed_reader = csv.reader(bed_file,delimiter='\t')

bed_writer_path = open(bed_path.split('.bed')[0]+'_Designed.with_NM.bed','w')
bed_writer = csv.writer(bed_writer_path,delimiter='\t')

###############
## METHODE 2 ##
###############

#custom_canonical_list = ['NM_000546']
amplnum = 1

bed_reader.next() # header
for bed_line in bed_reader:
	ampl_chr = bed_line[0]
	ampl_start = bed_line[1]
	ampl_stop = bed_line[2]
	try:
		amplicon = bed_line[3]
	except:
		amplicon = 'ampl' + str(amplnum)
		amplnum = amplnum + 1
	line2write = [ampl_chr,ampl_start,ampl_stop,amplicon,'0','+','.']
	NM = '???'
	gene = '???'
	exon = '??'
	ampl_pos = (int(bed_line[1])+int(bed_line[2]))/2 # position milieu ampl

	# FIND POSSIBLE CANONICAL USCS TRANSCRIPTS
	knownCanonical_file.seek(0) # begin after header
	knownCanonical_reader.next()
	for knownCanonical_line in knownCanonical_reader:
		canonical_start = int(knownCanonical_line[1])
		canonical_stop = int(knownCanonical_line[2])
		if (knownCanonical_line[0] == ampl_chr) and (canonical_start <= ampl_pos and canonical_stop >= ampl_pos) :
			canonical = knownCanonical_line[4]
			# FIND CORRESPONDING NM
			knownToRefSeq_file.seek(0) # begin after header
			knownToRefSeq_reader.next()
			for knownToRefSeq_line in knownToRefSeq_reader:
				if knownToRefSeq_line[0] == canonical:
					NM = knownToRefSeq_line[1]
					refGene_file.seek(0) # begin after header
					refGene_reader.next()
					for refgene_line in refGene_reader:
						if refgene_line[1] == NM :
							gene = refgene_line[12]
							strand = refgene_line[3]
							exonStarts = refgene_line[9].split(',')
							try:
								exonStarts.remove('')
							except:
								pass
							exonEnds = refgene_line[10].split(',')
							try:
								exonEnds.remove('')
							except:
								pass
							for index in range(len(exonStarts)):
								if (int(exonStarts[index]) < int(ampl_start) < int(exonEnds[index])) or (int(exonStarts[index]) < int(ampl_stop) < int(exonEnds[index])):
									if strand == '+':
										exon = index + 1
									if strand == '-':
										exon = len(exonStarts) - index
									break
							if exon == '??': #previous method didnt work, maybe overlapping amplicon
								for index in range(len(exonStarts)):
									if (int(ampl_start) < int(exonStarts[index])) and (int(ampl_stop) > int(exonEnds[index])):
										if strand == '+':
											exon = index + 1
										if strand == '-':
											exon = len(exonStarts) - index
										break
							break
					break
	gene_ex = gene + '_ex' + str(exon)
	line2write.append("GENE_ID=%s;SUBMITTED_REGION=%s" % (gene_ex,NM))
	bed_writer.writerow(line2write)
