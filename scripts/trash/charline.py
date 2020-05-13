#!/usr/bin/env python
import openpyxl
import glob
import csv

uv2pos = {}
unknow_variants = {}

unknow_xlsx = '/media/n06lbth/sauvegardes_pgm/SBT/Auto_user_S5-0198-82-NGS243_244_35pM_Chef_SBT-colon-lung_v7_530_194/IONXPRESS_029/new__IonXpress_029.finalReport.xlsx'
finalreport = openpyxl.load_workbook(unknow_xlsx)
annoSheet = finalreport.get_sheet_by_name('Annotation')
annoSheet_rows = tuple(annoSheet.rows)
for i in range(1,len(annoSheet_rows)):
	chrom = annoSheet_rows[i][1].value
	gene = annoSheet_rows[i][3].value
	c = annoSheet_rows[i][5].value
	pos = annoSheet_rows[i][15].value
	if chrom != None:
		unknow_variants[(chrom,gene,c)] = []
		uv2pos[(chrom,gene,c)] = pos


#xlsx_list = glob.glob('/media/n06lbth/sauvegardes_pgm/SBT/Auto_user_S5-0198-80-Run239_240_35pM_Chef_SBT-colon-lung_v7_530_192/*/*.xlsx')	
xlsx_list = glob.glob('/media/n06lbth/sauvegardes_pgm/SBT/Auto_user_S5-0198-78-Run237_238_35pM_Chef_SBT-colon-lung_v7_530_190/*/*.xlsx')
for xlsx in xlsx_list:
	name = xlsx.split('/')[-1].replace('.finalReport.xlsx','')
	finalreport = openpyxl.load_workbook(xlsx)
	try:
		annoSheet = finalreport.get_sheet_by_name('Annotation')
		annoSheet_rows = tuple(annoSheet.rows)
		for i in range(1,len(annoSheet_rows)):
			chrom = annoSheet_rows[i][1].value
			gene = annoSheet_rows[i][3].value
			c = annoSheet_rows[i][5].value
			thisvar = (chrom,gene,c)
			if thisvar in unknow_variants.keys():
				unknow_variants[thisvar].append(name)
	except:
		print 'dont work for %s' % xlsx

csvfile = open('/media/stuff/charline.tsv','w')
results = csv.writer(csvfile, delimiter ='\t')
for var in unknow_variants.keys():
	liste = [uv2pos[var],var[0],var[1],var[2]]
	for item in unknow_variants[var]:
		liste.append(item)
	results.writerow(liste)
