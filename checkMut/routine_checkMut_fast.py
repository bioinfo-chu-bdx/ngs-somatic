#!/usr/bin/env python
import os
import sys
import time
import subprocess

run_folder = sys.argv[1]
pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
checkMut_path = '%s/checkMut/checkMut.py' % pipeline_folder
FNULL = open(os.devnull, 'w')

# DIVERS MUTATIONS PAR DEFAUT
print "[%s] various SBT checkMut (9)..." % time.strftime("%H:%M:%S")
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','ERBB2','--cpos','c.2313_2324dup','--run-type','SBT'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','PIK3CA','--cpos','c.3140A>G','--run-type','SBT'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','PIK3CA','--cpos','c.3204_3205insA','--run-type','SBT'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','STK11','--cpos','c.169del','--run-type','SBT'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','STK11','--cpos','c.169dup','--run-type','SBT'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','STK11','--cpos','c.580G>T','--run-type','SBT'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','STK11','--cpos','c.595G>T','--run-type','SBT'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','FGFR3','--cpos','c.850del','--run-type','SBT'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','KIT','--cpos','c.2447A>T','--run-type','SBT'],stdout=FNULL,stderr=FNULL)

# LAM HOTSPOTS HOMOPOLYMERES
print "[%s] various LAM checkMut (12)..." % time.strftime("%H:%M:%S")
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','ASXL1','--cpos','c.2535dup','--run-type','LAM'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','ASXL1','--cpos','c.1927dup','--run-type','LAM'],stdout=FNULL,stderr=FNULL) # dupG ASXL1
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','ASXL1','--cpos','c.1927del','--run-type','LAM'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','ASXL1','--cpos','c.4127dup','--run-type','LAM'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','ASXL1','--cpos','c.4127del','--run-type','LAM'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','GATA2','--cpos','c.599del','--run-type','LAM'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','GATA2','--cpos','c.599dup','--run-type','LAM'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','GATA2','--cpos','c.302del','--run-type','LAM'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','SRSF2','--cpos','c.287del','--run-type','LAM'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','TP53','--cpos','c.455dup','--run-type','LAM'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','TP53','--cpos','c.455del','--run-type','LAM'],stdout=FNULL,stderr=FNULL)
subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','TP53','--cpos','c.454_466del','--run-type','LAM'],stdout=FNULL,stderr=FNULL)

# LISTE MUTATIONS ABL1
abl1_c_list = [
'c.944C>T',
'c.949T>C',
'c.951C>A',
'c.951C>G',
'c.667C>T',
'c.730A>G',
'c.758A>T',
'c.757T>C',
'c.763G>A',
'c.764A>T',
'c.880A>G',
'c.895G>C',
'c.943A>G',
'c.944_945delinsTG',
'c.949T>G',
'c.949T>A',
'c.950T>G',
'c.1010C>T',
'c.1075T>G',
'c.1076T>G',
'c.1075T>A',
'c.1392C>G',
'c.1393C>T',
'c.1402G>T',
'c.1504A>C'
]

print "[%s] ABL1 checkMut (%s)..." % (time.strftime("%H:%M:%S"),len(abl1_c_list))
for c_pos in abl1_c_list:
	subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','ABL1','--cpos',c_pos,'--run-type','ABL1','--sub-folder','ABL1'],stdout=FNULL,stderr=FNULL)
	
# LISTE MUTATIONS TABLEAU ADN CIRCULANT AUDREY
egfr_c_list = [
'c.2156G>C',
'c.2235_2249del',
'c.2236_2250del',
'c.2240_2257del',
'c.2236_2253del',
'c.2237_2255delinsT',
'c.2239_2256del',
'c.2239_2248delinsC',
'c.2240_2254del',
'c.2239_2256delinsCAA',
'c.2237_2252delinsT',
'c.2573T>G',
'c.2573_2574delinsGT',
'c.2582T>A',
'c.2369C>T',
'c.2389T>A',
'c.2390G>C',
'c.2303G>T',
'c.2307_2308insGCCAGCGTG',
'c.2309_2310delinsCCAGCGTGGAT',
'c.2310_2311insGGT',
'c.2319_2320insCAC',
'c.2311_2312insGCGTGGACA',
'c.2389T>A', # rajout samuel
'c.2390G>C'	# rajout samuel
]

print "[%s] EGFR checkMut (%s)..." % (time.strftime("%H:%M:%S"),len(egfr_c_list))
for c_pos in egfr_c_list:
	subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','EGFR','--cpos',c_pos,'--run-type','SBT','--sub-folder','EGFR'],stdout=FNULL,stderr=FNULL)
		
# LISTE MUTATIONS ADNc BRAF, NRAS, KRAS CHARLINE
braf_c_list = [
'c.1397G>C',
'c.1397G>T',
'c.1406G>C',
'c.1405G>A',
'c.1406G>T',
'c.1799T>A',
'c.1799_1800delinsAA',
'c.1799_1800delinsAT',
'c.1798_1799delinsAA',
'c.1798_1799delinsAG',
'c.1797_1799delinsGAG',
'c.1801A>G'
]

print "[%s] BRAF checkMut (%s)..." % (time.strftime("%H:%M:%S"),len(braf_c_list))
for c_pos in braf_c_list:
	subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','BRAF','--cpos',c_pos,'--run-type','SBT','--sub-folder','BRAF'],stdout=FNULL,stderr=FNULL)
	
nras_c_list = [
'c.35G>C',
'c.34G>T',
'c.35G>A',
'c.34G>C',
'c.34G>A',
'c.35G>T',
'c.38G>C',
'c.37G>T',
'c.38G>A',
'c.37G>C',
'c.37G>A',
'c.38G>T',
'c.52G>A',
'c.183A>T',
'c.183A>C',
'c.181C>A',
'c.182A>T',
'c.182A>C',
'c.182A>G',
'c.176C>A',
'c.175G>A',
'c.351G>C',
'c.351G>T',
'c.436G>A',
'c.437C>T'
]

print "[%s] NRAS checkMut (%s)..." % (time.strftime("%H:%M:%S"),len(nras_c_list))
for c_pos in nras_c_list:
	subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','NRAS','--cpos',c_pos,'--run-type','SBT','--sub-folder','NRAS'],stdout=FNULL,stderr=FNULL)
	
kras_c_list = [
'c.35G>C',
'c.34G>T',
'c.35G>A',
'c.34G>C',
'c.34G>A',
'c.35G>T',
'c.38G>C',
'c.37G>T',
'c.38G>A',
'c.37G>C',
'c.37G>A',
'c.38G>T',
'c.176C>A',
'c.176C>G',
'c.175G>T',
'c.175G>A',
'c.181C>G',
'c.183A>C',
'c.183A>T',
'c.181C>A',
'c.182A>T',
'c.182A>C',
'c.182A>G',
'c.351A>C',
'c.351A>T',
'c.436G>C',
'c.436G>A',
'c.437C>T',
]

print "[%s] KRAS checkMut (%s)..." % (time.strftime("%H:%M:%S"),len(kras_c_list))
for c_pos in kras_c_list:
	subprocess.Popen(['python',checkMut_path,'--run-folder',run_folder,'--gene','KRAS','--cpos',c_pos,'--run-type','SBT','--sub-folder','KRAS'],stdout=FNULL,stderr=FNULL)

pid_found = True
while pid_found:
	time.sleep(30)
	try:
		pid_found = subprocess.check_output(['pgrep','-f','%s/checkMut/checkMut.py' % pipeline_folder])
		print "[%s] \t running pid : %s" % (time.strftime("%H:%M:%S"),pid_found.replace('\n',','))
	except subprocess.CalledProcessError, e:
		pid_found = False
print "[%s] Done." % time.strftime("%H:%M:%S")
