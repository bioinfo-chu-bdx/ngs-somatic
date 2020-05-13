#!/usr/bin/python
import subprocess

run_to_clean = '/media/n06lbth/sauvegardes_pgm/SBT/Run_400-499/Auto_user_S5-0198-365-NGS461_462_35pM_Chef_SBT-colon-lung_v10_530_578/'

subprocess.Popen('rm %s/*/*processed*'%run_to_clean,shell=True)
subprocess.Popen('rm %s/*/*finalReport*'%run_to_clean,shell=True)
subprocess.Popen('rm %s/*/intermediate_files.zip'%run_to_clean,shell=True)
subprocess.Popen('rm -r %s/*/*intermediate_files'%run_to_clean,shell=True)
subprocess.Popen('rm %s/*/*.vbs'%run_to_clean,shell=True)
subprocess.Popen('rm -r %s/_plotCoverage'%run_to_clean,shell=True)
subprocess.Popen('rm -r %s/_dupG'%run_to_clean,shell=True)
subprocess.Popen('rm -r %s/_CNA'%run_to_clean,shell=True)
subprocess.Popen('rm -r %s/_checkMut'%run_to_clean,shell=True)
subprocess.Popen('rm -r %s/_Check-contamination'%run_to_clean,shell=True)
subprocess.Popen('rm %s/*.vbs'%run_to_clean,shell=True)
