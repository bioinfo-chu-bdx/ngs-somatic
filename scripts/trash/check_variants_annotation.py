#!/usr/bin/env python
import sys
import csv

alleles_file_path = sys.argv[1]
tsv_file_path = sys.argv[2]
	
alleles_file = open(alleles_file_path,"r")
alleles_file_reader = csv.reader(alleles_file,delimiter="\t")
tsv_file = open(tsv_file_path,"r")
tsv_file_reader = csv.reader(tsv_file,delimiter="\t")

missing_count = 0
variants_uniq = []

alleles_file_reader.next()
for alleles_line in alleles_file_reader:
	if alleles_line[4] == "Homozygous" or alleles_line[4] == "Heterozygous":
		variant = (alleles_line[1],alleles_line[2],alleles_line[3])# pos, ref, alt
		if variant not in variants_uniq:
			variants_uniq.append(variant)
			found = False
			tsv_file.seek(0)
			tsv_file_reader.next()
			for tsv_line in tsv_file_reader:
				if variant == (tsv_line[13],tsv_line[14],tsv_line[15]):
					print " - OK - " + str(variant)
					found = True
					continue
				elif variant == str(int(tsv_line[13])+1)+","+tsv_line[14]+","+tsv_line[15]:
					print " - OK - " + str(variant) + "  *warning* position shift -1"
					found = True
					continue
			if not found:
				print " - NOK!!! - " + str(variant)
				missing_count = missing_count+1
		else:
			print "*warning* ignoring duplicate " + variant
alleles_file.close()
tsv_file.close()
print "\n######## SUMMARY ########\n --> MISSING VARIANTS IN VARIANTANNOTATION : " + str(missing_count)
