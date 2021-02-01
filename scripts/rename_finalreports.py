#!/usr/bin/python
import os
import sys
import glob
import subprocess

run = '/media/n06lbth/sauvegardes_pgm/SBT/Auto_user_S5-0198-505-NGS538_539_35pM_Chef_SBT-colon-lung_v10_530_736'

fps = glob.glob('%s/*/*.finalReport.xlsx' % run)
for fp in fps:
	renamed_fp = fp.replace('finalReport.xlsx','finalReport.new.xlsx')
	subprocess.call(['mv',fp,renamed_fp])
