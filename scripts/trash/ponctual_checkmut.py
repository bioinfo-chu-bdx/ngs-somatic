#!/usr/bin/env python
import sys
import subprocess

checkMut_path = '/DATA/work/scripts/mpileup_checkMut.py'

run_folders = [
	'/media/n06lbth/sauvegardes_pgm/SBT/Run_400-499/Auto_user_S5-0198-302-NGS410_411_35pM_Chef_SBT-colon-lung_v10_530_509/',
	'/media/n06lbth/sauvegardes_pgm/SBT/Run_400-499/Auto_user_S5-0198-305-NGS412_413_35p_Lib21pB_30p_Chef_SBT-colon-lung_v10_530_514/',
	'/media/n06lbth/sauvegardes_pgm/SBT/Run_400-499/Auto_user_S5-0198-307-NGS414_415_35pM_Chef_SBT-colon-lung_v10_530_516/',
	'/media/n06lbth/sauvegardes_pgm/SBT/Run_400-499/Auto_user_S5-0198-309-Run416_417_35pM_Run_lib22_B_30pM_Chef_SBT-colon-lung_v10_530_520/',
	'/media/n06lbth/sauvegardes_pgm/SBT/Auto_user_S5-0198-311-NGS418_419_35pM_Chef_SBT-colon-lung_v10_530_523/',
]

kras_c_list = [
	'c.34G>T ',
	'c.34G>C ',
	'c.34G>A ',
	'c.35G>T ',
	'c.35G>A ',
	'c.35G>C',
	'c.38G>A ',
]

for run_folder in run_folders:
	for c_pos in kras_c_list:
		subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','KRAS','--cpos',c_pos,'--run-type','SBT','--sub-folder','KRAS_G12_G13_for_NTC'])
