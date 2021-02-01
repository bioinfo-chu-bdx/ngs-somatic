#!/usr/bin/env python
import os
import sys
import time
import json
import sqlite3
import subprocess


def dict_factory(cursor, row):
	d = {}
	for idx, col in enumerate(cursor.description):
		d[col[0]] = row[idx]
	return d

run_folder = sys.argv[1]
pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
checkMut_path = '%s/checkMut/checkMut.py' % pipeline_folder
FNULL = open(os.devnull, 'w')

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.loads(g.read().replace('$NGS_PIPELINE_BX_DIR',os.environ['NGS_PIPELINE_BX_DIR']))

db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

if run_folder.endswith('/'):
	run_folder = run_folder[:-1]
run_name = os.path.basename(run_folder)

db_cur.execute("SELECT platform FROM Run WHERE runID='%s'" % run_name)
db_run = db_cur.fetchone()
platform = db_run['platform']
print "- running routine checkMuts for %s platform" % platform

if platform.lower() == 'illumina':
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.23A>C','--run-type','LAM']) # '--min-sample','30'
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.29A>C','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.292A>C','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.323A>C','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.722T>G','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','RUNX1','--cpos','c.1265A>C','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','RUNX1','--cpos','c.1270T>G','--run-type','LAM'])
	# subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.564_566del','--run-type','LAM'])
	# subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.688A>C','--run-type','LAM'])
	# subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.700A>T','--run-type','LAM'])
	# subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.737T>G','--run-type','LAM'])
	# subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.747T>G','--run-type','LAM'])
	# subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.751A>C','--run-type','LAM'])
	# subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','KMT2C','--cpos','c.2656C>T','--run-type','LAM'])
	# subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','KMT2C','--cpos','c.5053G>T','--run-type','LAM'])
	# subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.230_233dup','--run-type','LAM'])
	# subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','CEBPA','--cpos','c.992T>A','--run-type','LAM'])

elif platform.lower() == 'ion torrent':

	# DIVERS MUTATIONS PAR DEFAUT
	print "[%s] various SBT checkMut (9)..." % time.strftime("%H:%M:%S")
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','ERBB2','--cpos','c.2313_2324dup','--run-type','SBT'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','PIK3CA','--cpos','c.3140A>G','--run-type','SBT'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','PIK3CA','--cpos','c.3204_3205insA','--run-type','SBT'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','STK11','--cpos','c.169del','--run-type','SBT'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','STK11','--cpos','c.169dup','--run-type','SBT'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','STK11','--cpos','c.580G>T','--run-type','SBT'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','STK11','--cpos','c.595G>T','--run-type','SBT'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','FGFR3','--cpos','c.850del','--run-type','SBT'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','KIT','--cpos','c.2447A>T','--run-type','SBT'])

	# LAM HOTSPOTS HOMOPOLYMERES
	print "[%s] various LAM checkMut (12)..." % time.strftime("%H:%M:%S")
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','ASXL1','--cpos','c.2535dup','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','ASXL1','--cpos','c.1927dup','--run-type','LAM']) # dupG ASXL1
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','ASXL1','--cpos','c.1927del','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','ASXL1','--cpos','c.4127dup','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','ASXL1','--cpos','c.4127del','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','GATA2','--cpos','c.599del','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','GATA2','--cpos','c.599dup','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','GATA2','--cpos','c.302del','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','SRSF2','--cpos','c.287del','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','TP53','--cpos','c.455dup','--run-type','TP53,LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','TP53','--cpos','c.455del','--run-type','TP53,LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','TP53','--cpos','c.454_466del','--run-type','TP53,LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','TP53','--cpos','c.743G>A','--run-type','TP53,LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','TP53','--cpos','c.796G>A','--run-type','TP53,LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','TP53','--cpos','c.817C>T','--run-type','TP53,LAM'])
	## ILLUMINA
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','CEBPA','--cpos','c.23A>C','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','CEBPA','--cpos','c.29A>C','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','CEBPA','--cpos','c.292A>C','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','CEBPA','--cpos','c.323A>C','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','CEBPA','--cpos','c.564_566del','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','CEBPA','--cpos','c.688A>C','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','CEBPA','--cpos','c.700A>T','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','CEBPA','--cpos','c.722T>G','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','CEBPA','--cpos','c.737T>G','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','CEBPA','--cpos','c.747T>G','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','CEBPA','--cpos','c.751A>C','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','RUNX1','--cpos','c.1265A>C','--run-type','LAM'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','RUNX1','--cpos','c.1270T>G','--run-type','LAM'])

	# DIVERS LYMPHOME B
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','MYD88','--cpos','c.794T>C','--run-type','Lymphome_B'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','NOTCH2','--cpos','c.7198C>T','--run-type','Lymphome_B'])
	subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','XPO1','--cpos','c.1711G>A','--run-type','Lymphome_B'])

	# LISTE MUTATIONS ABL1
	abl1_c_list = [
	'c.944C>T',
	'c.949T>C',
	'c.951C>A',
	'c.951C>G',
	'c.667C>T',
	'c.730A>G',
	'c.749G>A', # ajout lisa 15/05/2020
	'c.758A>T',
	'c.757T>C',
	'c.763G>A',
	'c.764A>T',
	'c.880A>G',
	'c.895G>C',
	'c.895G>T',
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
		subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--gene','ABL1','--cpos',c_pos,'--run-type','ABL1'])#,'--sub-folder','ABL1'])
		
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
		subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','EGFR','--cpos',c_pos,'--run-type','SBT'])#,'--sub-folder','EGFR'])
			
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
	'c.1801A>G',
	'c.1797delinsTACTACG', # ajout cq myel
	'c.1798_1799delinsCA',
	'c.1798G>A',
	'c.1798G>C',
	'c.1798G>T',
	'c.1799_1800delinsAA',
	'c.1799_1801del',
	'c.1799T>C',
	'c.1799T>G',
	'c.1801_1803del',
	'c.1802A>T',
	'c.1803A>C',
	'c.1803A>T'
	]

	print "[%s] BRAF checkMut (%s)..." % (time.strftime("%H:%M:%S"),len(braf_c_list))
	for c_pos in braf_c_list:
		subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','BRAF','--cpos',c_pos,'--run-type','SBT,Lymphome_B'])#,'--sub-folder','BRAF'])
		
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
		subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','NRAS','--cpos',c_pos,'--run-type','SBT'])#,'--sub-folder','NRAS'])
		
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
		subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','KRAS','--cpos',c_pos,'--run-type','SBT'])#,'--sub-folder','KRAS'])
		
	cd79b_c_list = [
	'c.591G>C',
	'c.589G>A',
	'c.588C>T',
	'c.588C>G',
	'c.587A>T',
	'c.587A>G',
	'c.587A>C',
	'c.586T>G',
	'c.586T>C',
	'c.586T>A'
	]

	print "[%s] CD79B checkMut (%s)..." % (time.strftime("%H:%M:%S"),len(kras_c_list))
	for c_pos in cd79b_c_list:
		subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','CD79B','--cpos',c_pos,'--run-type','Lymphome_B'])#,'--sub-folder','KRAS'])

	stat3_c_list = [
	'c.1771A>G',
	'c.1772A>T',
	'c.1778G>C',
	'c.1825A>G',
	'c.1831A>G',
	'c.1832G>A',
	'c.1832G>T',
	'c.18040A>G',
	'c.1842C>G',
	'c.1847A>G',
	'c.1846G>A',
	'c.1850G>A',
	'c.1850G>T',
	'c.1849G>A',
	'c.1852G>C',
	'c.1853G>A',
	'c.1858A>G',
	'c.1859C>G',
	'c.1861T>G',
	'c.1863C>G',
	'c.1862T>C',
	'c.1865C>T',
	'c.1881C>A',
	'c.1907C>T',
	'c.1907C>A',
	'c.1909G>A',
	'c.1909G>T',
	'c.1910T>C',
	'c.1913A>G',
	'c.1915C>G',
	'c.1915C>T',
	'c.1915C>A',
	'c.1919A>T',
	'c.1924A>G',
	'c.1927C>A',
	'c.1929A>C',
	'c.1931_1933del',
	'c.1938C>G',
	'c.1939A>G',
	'c.1940A>T',
	'c.1954G>A',
	'c.1967G>A',
	'c.1968delinsTTTT',
	'c.1970A>G',
	'c.1970A>C',
	'c.1969T>A',
	'c.1969_1971dup',
	'c.1969_1980dup',
	'c.1973A>T',
	'c.1974G>T',
	'c.1972_1974delinsTAT',
	'c.1972A>G',
	'c.1976A>T',
	'c.1975A>C',
	'c.1979T>G',
	'c.1978T>A',
	'c.1981_1982delinsAT',
	'c.1981G>T',
	'c.1981A>T',
	'c.1981G>C',
	'c.1985C>T',
	'c.1988_1989delinsTT',
	'c.1998T>A',
	'c.1999C>G',
	'c.2003C>T',
	'c.2003C>A',
	'c.2104G>A'
	]

	print "[%s] STAT3 checkMut (%s)..." % (time.strftime("%H:%M:%S"),len(stat3_c_list))
	for c_pos in stat3_c_list:
		subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','STAT3','--cpos',c_pos,'--run-type','Lymphome_T','--sub-folder','STAT3'])#,'--sub-folder','KRAS'])

	stat5b_c_list = [
	'c.1787G>T',
	'c.1883C>G',
	'c.1888G>C',
	'c.1901A>T',
	'c.1907A>C',
	'c.1924A>C',
	'c.1937T>C',
	'c.1942A>T',
	'c.1975C>T',
	'c.1994A>T',
	'c.1993T>C',
	'c.2110A>C',
	'c.1735G>A'
	]

	print "[%s] STAT5B checkMut (%s)..." % (time.strftime("%H:%M:%S"),len(stat5b_c_list))
	for c_pos in stat5b_c_list:
		subprocess.call(['python',checkMut_path,'--run-folder',run_folder,'--min-sample','30','--gene','STAT5B','--cpos',c_pos,'--run-type','Lymphome_T','--sub-folder','STAT5B'])#,'--sub-folder','KRAS'])

print "[%s] Done." % time.strftime("%H:%M:%S")
