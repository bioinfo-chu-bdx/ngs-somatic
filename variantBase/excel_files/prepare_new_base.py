#!/usr/bin/env python
# coding: utf-8
import openpyxl
import sys

variantBase_empty_path = sys.argv[1]
#gene_list = sys.argv[2]
if len(sys.argv) == 3:
	bed_with_nm = sys.argv[2]
	list_bed = [bed_with_nm]
else:
	list_bed = sys.argv[2:]

glist = []

for bed_with_nm in list_bed:
	target_bed = open(bed_with_nm,'r')
	target_bed.next()
	for line in target_bed:
		gene = line.split('GENE_ID=')[-1].split(';')[0].split('_')[0]
		if gene not in glist :
			glist.append(gene)

vb_empty = openpyxl.load_workbook(variantBase_empty_path)
gene_sheet = vb_empty['Gene']

for gene in glist:
	newsheet = vb_empty.copy_worksheet(gene_sheet)
	newsheet.title = gene

vb_empty.remove_sheet(gene_sheet)
vb_empty.save('/DATA/work/variantBase/newBase_files/VariantBase_%s.xlsx' % sys.argv[2].split('/')[-1])
