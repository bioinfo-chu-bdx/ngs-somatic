#!/usr/bin/env python
import os
import sys
import time
import subprocess

run_folder = sys.argv[1]
pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
checkMut_path = '%s/checkMut/checkMut.py' % pipeline_folder
FNULL = open(os.devnull, 'w')

## ILLUMINA
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.23A>C','--run-type','LAM']) # '--min-sample','30'
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.29A>C','--run-type','LAM'])
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.292A>C','--run-type','LAM'])
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.323A>C','--run-type','LAM'])
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.564_566del','--run-type','LAM'])
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.688A>C','--run-type','LAM'])
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.700A>T','--run-type','LAM'])
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.722T>G','--run-type','LAM'])
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.737T>G','--run-type','LAM'])
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.747T>G','--run-type','LAM'])
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.751A>C','--run-type','LAM'])
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','RUNX1','--cpos','c.1265A>C','--run-type','LAM'])
subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','RUNX1','--cpos','c.1270T>G','--run-type','LAM'])

print "[%s] Done." % time.strftime("%H:%M:%S")
