#!/usr/bin/env python
import os
import sys
import openpyxl
import csv

temoincnv_paths = [
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-182-Run_lib11_50pM_Lymphome_T_530_351/temoinCNV1-EFS396/temoinCNV1-EFS396_IonXpress_030.finalReport.xlsx',
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-182-Run_lib11_50pM_Lymphome_T_530_351/temoinCNV2-EFS218/temoinCNV2-EFS218_IonXpress_031.finalReport.xlsx',
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-182-Run_lib11_50pM_Lymphome_T_530_351/temoinCNV3-EFS307/temoinCNV3-EFS307_IonXpress_032.finalReport.xlsx',
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-183-Run_lib12_50pM_Lymphome_T_530_352/temoinCNV4-EFS362/temoinCNV4-EFS362_IonXpress_034.finalReport.xlsx',
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-183-Run_lib12_50pM_Lymphome_T_530_352/temoinCNV5-EFS293/temoinCNV5-EFS293_IonXpress_035.finalReport.xlsx',
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-184-Run_lib13_50pM_Lymphome_T_530_353/temoinCNV6-EFS331/temoinCNV6-EFS331_IonXpress_045.finalReport.xlsx',
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-184-Run_lib13_50pM_Lymphome_T_530_353/temoinCNV7-EFST109h/temoinCNV7-EFST109h_IonXpress_046.finalReport.xlsx'
]

variants_seen = {}
pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
csvresultsfile = open('%s/scripts/temoinCNV_variants_lymphomeT.tsv' % pipeline_folder,'w')
results = csv.writer(csvresultsfile, delimiter = '\t')

for temoincnv_path in temoincnv_paths:
	control_name = temoincnv_path.split('/')[-1].split('_IonXpress')[0]
	temoincnv_report = openpyxl.load_workbook(temoincnv_path)
	annoSheet = temoincnv_report.get_sheet_by_name('Annotation')
	annoSheet_rows = tuple(annoSheet.rows)
	for i in range(1,len(annoSheet_rows)):
		nm = annoSheet_rows[i][2].value
		c = annoSheet_rows[i][5].value
		if not nm:
			continue
		if (nm,c) in variants_seen:
			variants_seen[(nm,c)] = variants_seen[(nm,c)] + ',' + control_name
		else:
			variants_seen[(nm,c)] = control_name
			
row_num = 0			
for variant in variants_seen:
	results.writerow([variant[0],variant[1],variants_seen[variant]])
		
csvresultsfile.close()
