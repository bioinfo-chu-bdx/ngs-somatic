#!/usr/bin/env python
import glob
import json
import sys
import os
import subprocess

run_folder = sys.argv[1]
cna_dir = run_folder + '/_CNA'
fnull = open(os.devnull, 'w')

oncocontrols_B = [
'/media/n06lbth/sauvegardes_pgm/Lymphome_B/Auto_user_S5-0198-197-lib8_aout2018_30pM_Lymphome_B_Oceane_SBT_530_372/temoinCNV1-EFS396/temoinCNV1-EFS396_IonXpress_086.bam',
'/media/n06lbth/sauvegardes_pgm/Lymphome_B/Auto_user_S5-0198-197-lib8_aout2018_30pM_Lymphome_B_Oceane_SBT_530_372/temoinCNV2-EFS218/temoinCNV2-EFS218_IonXpress_087.bam',
'/media/n06lbth/sauvegardes_pgm/Lymphome_B/Auto_user_S5-0198-197-lib8_aout2018_30pM_Lymphome_B_Oceane_SBT_530_372/temoinCNV3-EFS307/temoinCNV3-EFS307_IonXpress_088.bam',
'/media/n06lbth/sauvegardes_pgm/Lymphome_B/Auto_user_S5-0198-197-lib8_aout2018_30pM_Lymphome_B_Oceane_SBT_530_372/temoinCNV4-EFS362/temoinCNV4-EFS362_IonXpress_089.bam',
'/media/n06lbth/sauvegardes_pgm/Lymphome_B/Auto_user_S5-0198-198-lib9_aout2018_30pM_Lymphome_B_Oceane_SBT_530_373/temoinCNV5-EFS293/temoinCNV5-EFS293_IonXpress_002.bam',
'/media/n06lbth/sauvegardes_pgm/Lymphome_B/Auto_user_S5-0198-200-lib10_aout2018_30pM_Lymphome_B_Oceane_SBT_530_374/temoinCNV6-EFS331/temoinCNV6-EFS331_IonXpress_007.bam'
]
oncocontrols_T = [
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-182-Run_lib11_50pM_Lymphome_T_530_351/temoinCNV1-EFS396/temoinCNV1-EFS396_IonXpress_030.bam',
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-182-Run_lib11_50pM_Lymphome_T_530_351/temoinCNV2-EFS218/temoinCNV2-EFS218_IonXpress_031.bam',
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-182-Run_lib11_50pM_Lymphome_T_530_351/temoinCNV3-EFS307/temoinCNV3-EFS307_IonXpress_032.bam',
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-183-Run_lib12_50pM_Lymphome_T_530_352/temoinCNV4-EFS362/temoinCNV4-EFS362_IonXpress_034.bam',
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-183-Run_lib12_50pM_Lymphome_T_530_352/temoinCNV5-EFS293/temoinCNV5-EFS293_IonXpress_035.bam',
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-184-Run_lib13_50pM_Lymphome_T_530_353/temoinCNV6-EFS331/temoinCNV6-EFS331_IonXpress_045.bam',
'/media/n06lbth/sauvegardes_pgm/Lymphome_T/Auto_user_S5-0198-184-Run_lib13_50pM_Lymphome_T_530_353/temoinCNV7-EFST109h/temoinCNV7-EFST109h_IonXpress_046.bam'
]

# RECUPERATION BAM ET RUNTYPES 
bamlist = glob.glob(run_folder+'/*/*.bam')
bamlist = [item for item in bamlist if not 'processed' in item]

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

if os.path.isfile(run_folder+'/barcodes.json'):
	with open(run_folder+'/barcodes.json', 'r') as g:
		barcodes_json = json.load(g)

runtypes = []
oncotests_B = []
oncotests_T = []
for bamfile in bamlist:
	sample = bamfile.split('/')[-1].split('_IonXpress')[0]
	barcode = 'IonXpress_' + bamfile.split('IonXpress_')[-1].split('.bam')[0]
	target = barcodes_json[barcode]['target_region_filepath'].split('/')[-1]
	for runtype in global_param['run_type']:
		if global_param['run_type'][runtype]['target_bed'].split('/')[-1] == target:
			if runtype not in runtypes:
				runtypes.append(runtype)
			break
	iscontrol = False
	for control_name in ['NTC_','H2O','ACROMETRIX','HD802','TEMOINCNV']:
		if control_name in sample.upper():
			iscontrol = True
			break
	if runtype == 'Lymphome_B' and not iscontrol:
		oncotests_B.append(bamfile)
	if runtype == 'Lymphome_T' and not iscontrol:
		oncotests_T.append(bamfile)

# CREATION DOSSIER _CNA
if not os.path.isdir(cna_dir):
	subprocess.call(['mkdir',cna_dir])

# CREATION DOSSIER _CNA/RUNTYPE
for runtype in runtypes:
	cna_runtype_dir = '%s/%s' % (cna_dir,runtype)
	if not os.path.isdir(cna_runtype_dir):
		subprocess.call(['mkdir',cna_runtype_dir])

####################
### CNV ANALYSIS ###
####################

subprocess.call(['python','%s/CNV/ioncopy/make_ioncopy_input.py' % pipeline_folder,'--run-folder', run_folder], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

for runtype in runtypes:
	cna_runtype_dir = '%s/%s/' % (cna_dir,runtype)
	if not os.path.isdir(cna_runtype_dir):
		subprocess.call(['mkdir',cna_runtype_dir])
		
	print "\t - ioncopy... (%s)" % runtype
	subprocess.call(['Rscript','--vanilla','%s/CNV/ioncopy/run_ioncopy.R' % pipeline_folder,cna_runtype_dir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	
	if runtype == 'Lymphome_B':
		print "\t - oncocnv... (%s)" % runtype
		oncobed = '%s/CNV/oncocnv/Lymphome_B_IAD119887_231.oncocnv.bed' % pipeline_folder
		subprocess.call(['mkdir',cna_runtype_dir+'/oncocnv'])
		controllist = ','.join(oncocontrols_B)
		testlist = ','.join(oncotests_B)
		subprocess.call(['bash','%s/CNV/oncocnv/ONCOCNV.sh' % pipeline_folder,controllist,testlist,oncobed,cna_runtype_dir+'/oncocnv'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		subprocess.call(['python','%s/CNV/make_cnv_finalreport_with_oncocnv.py' % pipeline_folder,cna_runtype_dir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
	elif runtype == 'Lymphome_T':
		print "\t - oncocnv... (%s)" % runtype
		oncobed = '%s/CNV/oncocnv/Lymphome_T_IAD120574_238.oncocnv.bed' % pipeline_folder
		subprocess.call(['mkdir',cna_runtype_dir+'/oncocnv'])
		controllist = ','.join(oncocontrols_T)
		testlist = ','.join(oncotests_T)
		subprocess.call(['bash','%s/CNV/oncocnv/ONCOCNV.sh' % pipeline_folder,controllist,testlist,oncobed,cna_runtype_dir+'/oncocnv'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		subprocess.call(['python','%s/CNV/make_cnv_finalreport_with_oncocnv.py' % pipeline_folder,cna_runtype_dir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
	else:
		subprocess.call(['python','%s/CNV/make_cnv_finalreport.py' % pipeline_folder,cna_runtype_dir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
