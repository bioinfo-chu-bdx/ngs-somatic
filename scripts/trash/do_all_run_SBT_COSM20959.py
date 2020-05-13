#!/usr/bin/env python
import subprocess
import time

list = open('/DATA/work/variantBase/temp/run_list_LYMPHOME_B.fullpath.txt','r')
for run in list:
	subprocess.call(['python','/DATA/work/scripts/all_COSM4781251.py','--run-folder',run.replace('\n','')])
