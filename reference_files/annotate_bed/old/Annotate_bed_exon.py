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

bed_path = sys.argv[1]#'/media/stuff/NM-Lisa.bed'
refGene_path = '/DATA/work/reference_files/Annotate_bed_script/refGene_sans_doublons_NM.txt'

refGene_file = open(refGene_path,'r')
bed_file = open(bed_path,'r')

refGene_reader = csv.reader(refGene_file,delimiter='\t')
bed_reader = csv.reader(bed_file,delimiter='\t')

bed_writer_path = open('/DATA/work/reference_files/Annotate_bed_script/bed_result.bed','w')
bed_writer = csv.writer(bed_writer_path,delimiter='\t')

###############
## METHODE 2 ##
###############

header = bed_reader.next()
bed_writer.writerow(header)
for bed_line in bed_reader:
	ampl_chr = bed_line[0]
	ampl_start = bed_line[1]
	ampl_stop = bed_line[2]
	amplicon = bed_line[3]
	NM = bed_line[7].split(';SUBMITTED_REGION=')[-1]
	gene = bed_line[7].split(';SUBMITTED_REGION=')[0].split('GENE_ID=')[-1].split('_ex')[0].split('_in')[0]
	line2write = [ampl_chr,ampl_start,ampl_stop,amplicon,'0','+','.']
	exon = '??'
	ampl_pos = (int(bed_line[1])+int(bed_line[2]))/2 # position milieu ampl

	refGene_file.seek(0) # begin after header
	refGene_reader.next()
	for refgene_line in refGene_reader:
		if refgene_line[1] == NM :
			geneee = refgene_line[12]
			if geneee != gene:
				print "PRBLM : LES NOMS DE GENE NE CORRESPONDENT PAS"
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
			
	gene_ex = gene + '_ex' + str(exon)
	line2write.append('GENE_ID=%s;SUBMITTED_REGION=%s' % (gene_ex,NM))
	bed_writer.writerow(line2write)
