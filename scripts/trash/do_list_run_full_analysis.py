#!/usr/bin/env python
import subprocess
import time

list = open('/DATA/work/variantBase/temp/run_list_LYMPHOME_T.fullpath.txt','r')
for run in list:
	subprocess.call(['python','/DATA/work/run_analysis.py','--full-run',run.replace('\n','')])
