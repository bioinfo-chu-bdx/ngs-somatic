#!/usr/bin/env python
import sys
import subprocess

run_folder_list = [
	'/media/n06lbth/sauvegardes_pgm/SBT/Run_100-199/Auto_user_PGM-250-Run159_CQ_35pM_Chef_SBT_colon_lung_v5_318v2_322/',
	'/media/n06lbth/sauvegardes_pgm/SBT/Run_100-199/Auto_user_S5-0198-48-NGS196_197_35pM_Chef_SBT-colon-lung_v7_530_156/',
	'/media/n06lbth/sauvegardes_pgm/SBT/Run_100-199/Auto_user_S5-0198-3-Run_validation_PGM159CQ160_40pM_Chef_SBT-colon-lung_v5_530_99/',
	'/media/n06lbth/sauvegardes_pgm/SBT/Run_200-299/Auto_user_S5-0198-95-Run255_257_35pM_Chef_SBT-colon-lung_v7_530_208/'
]
checkMut_path = '/DATA/work/scripts/mpileup_checkMut.py'

kras_c_list = [
'c.35G>C',
'c.34G>T',
'c.35G>A',
'c.34G>C',
'c.34G>A',
'c.35G>T',
'c.181C>G',
'c.183A>C',
'c.183A>T',
'c.181C>A',
'c.182A>T',
'c.182A>C',
'c.182A>G',
]

for run_folder in run_folder_list:
	for c_pos in kras_c_list:
		subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','KRAS','--cpos',c_pos])
