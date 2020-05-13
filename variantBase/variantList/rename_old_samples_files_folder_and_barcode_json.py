#!/usr/bin/python
import sys
import csv
import sqlite3
import zipfile
import json
import re
import os
import glob
import time
import uuid
from datetime import date
import subprocess

month2num = {'Jan':1,'Feb':2,'Mar':3,'Apr':4,'May':5,'Jun':6,'Jul':7,'Aug':8,'Sep':9,'Oct':10,'Nov':11,'Dec':12}
control_names = ['H2O','H20','NTC','ACROMETRIX','BAF5','BAF-5','HD300','HD301','HD802','HD901','HD748','HORIZON','TEMOIN'] # liste des noms possibles pour les temoins negatifs

run_list_path = '/DATA/work/variantBase/runList/runList_ALL.fullpath.txt'

#     __                   
#    |__) |  | |\ |  
#    |  \ \__/ | \|    	

sbt_pattern = r"^.*[A-Z]{2}[0-9]{3}[A-Z]{1}$" ## ex BESCOND-AW210F (2 lettres 3 chiffres 1 lettre)
sbt_old_pattern = r"^.*[A-Z]{1}[0-9]{3}[A-Z]{1}$" ## ex GUSMINI-Z492F (1 lettre 3 chiffres 1 lettre)
hemato_pattern = r"^.*[0-9]{3}\.[0-9]{3}$" ## ex ETIENNE-MICHELE-204.066 (se termine par 3 chiffres 1 point 3 chiffres)
hemato_old_pattern = r"^.*[0-9]{3}\.[0-9]{2}$" ## ex ETIENNE-MICHELE-204.06 (se termine par 3 chiffres 1 point 2 chiffres)
hemato_retard_pattern = r"^.*[0-9]{3}-[0-9]{3}$" ## ex ETIENNE-MICHELE-204-066 (se termine par 3 chiffres 1 point 3 chiffres) because fuck you
hemato_retard_old_pattern = r"^.*[0-9]{3}-[0-9]{2}$" ## ex ETIENNE-MICHELE-204-06 (se termine par 3 chiffres 1 point 2 chiffres) because fuck you
hemato_long_pattern = r"^.*TU\d\d\.c\d\.R\d\.B[0-9]{3}\.[0-9]{3}$" ## ex LACARRA-FRANCOIS-TU00.c0.R0.B201.069 (se termine par 3 chiffres 1 point 3 chiffres)

try:
	run_list = [sys.argv[1]]
except:
	run_list = []
	rl = open(run_list_path,'r')
	for run in rl:
		run_folder = run.replace('\n','')
		run_list.append(run_folder)
		
for run_folder in run_list:
	if run_folder.endswith('/'):
		run_name = run_folder.split('/')[-2]
	else:
		run_name = run_folder.split('/')[-1]
	print "#############################################################"
	print "RUN : %s" % run_name
	print "#############################################################"
	modifications = False
	with open(run_folder+'/barcodes.json', 'r') as g:
		barcodes_json = json.load(g)
		for barcode in barcodes_json:
			og_sample_name = barcodes_json[barcode]['sample']
			if 'sample_id' not in barcodes_json[barcode]:
				barcodes_json[barcode]['sample_id'] = ''
			og_sample_id = barcodes_json[barcode]['sample_id']
			sample_name = barcodes_json[barcode]['sample']

			# REPLACE SPACE WITH -, ELIMINATE END - or _
			sample_name = sample_name.replace(' ','-')
			sample_name = sample_name.replace('_','-')
			while sample_name.endswith('-') or sample_name.endswith('_'):
				sample_name = sample_name[:-1]
				
			#### REGEX pour detecter format SBT, ou HEMATO, ou aucun des deux et dans ce cas random id et le reste dans lastname
			print "-------------------------"
			print "* name -> %s" % (og_sample_name)
			
			id_in_filename = True
			iscontrol = False
			for control_name in control_names:
				if control_name in sample_name.upper():
					if re.match(sbt_pattern,sample_name) or re.match(sbt_old_pattern,sample_name):
						# control should not have dna number... this is to avoid REVLEG-AH200F
						continue
					else:
						iscontrol = True
			if iscontrol:
				print "\t -is control"
				random_uuid = uuid.uuid1()
				dna_number = 'CONTROL-%s' % random_uuid.hex[:8].upper()
				id_in_filename = False
			if not iscontrol:
				if re.match(sbt_pattern,sample_name):
					print "\t -sbt pattern match"
					dna_number = sample_name[-6:]
					sample_name = sample_name[:-7]
				elif re.match(sbt_old_pattern,sample_name):
					print "\t -sbt old pattern match"
					dna_number = sample_name[-5:]
					sample_name = sample_name[:-6]
				elif re.match(hemato_long_pattern,sample_name):
					print "\t -hemato long pattern match"
					dna_number = sample_name[-7:]
					sample_name = sample_name[:-20]
				elif re.match(hemato_pattern,sample_name):
					print "\t -hemato short pattern match"
					dna_number = sample_name[-7:]
					sample_name = sample_name[:-8]
				elif re.match(hemato_old_pattern,sample_name):
					print "\t -hemato short old pattern match"
					dna_number = sample_name[-6:]
					sample_name = sample_name[:-7]
				elif re.match(hemato_retard_pattern,sample_name):
					print "\t -hemato retarded pattern match"
					dna_number = sample_name[-7:].replace('-','.')
					sample_name = sample_name[:-8]
				elif re.match(hemato_retard_old_pattern,sample_name):
					print "\t -hemato retarded old pattern match"
					dna_number = sample_name[-6:].replace('-','.')
					sample_name = sample_name[:-7]
				else: # NO REGEX MATCH
					print "\t -no regex match"
					if barcodes_json[barcode]['sample_id'] == '':
						random_uuid = uuid.uuid1()
						random_uuid = 'ID'+random_uuid.hex[:8].upper()
						dna_number = random_uuid
					else:
						dna_number = barcodes_json[barcode]['sample_id']
					id_in_filename = False
				
			# IF DNA NUMBER = 000.000 -< unknow; generate uuid instead
			if dna_number == '000.000':
				if barcodes_json[barcode]['sample_id'] == '':
					random_uuid = uuid.uuid1()
					random_uuid = 'ID'+random_uuid.hex[:8].upper()
					dna_number = random_uuid
				else:
					dna_number = barcodes_json[barcode]['sample_id']
				id_in_filename = False
			
			# ELIMINATE sample_name END - or _
			while sample_name.endswith('-') or sample_name.endswith('_'):
				sample_name = sample_name[:-1]
			
			if id_in_filename:
				if sample_name != '':
					full_name = sample_name + '-' + dna_number
				else:
					full_name = dna_number
			else:
				full_name = sample_name
			if og_sample_id != dna_number:
				#ui = raw_input("* NEW SAMPLE ID = %s (OK?)" % dna_number) # MODIFY IN BARCODES JSON
				#if ui == '':
				print "* NEW SAMPLE ID = %s" % dna_number
				barcodes_json[barcode]['sample_id'] = dna_number
				modifications = True
			if og_sample_name != full_name:
				#ui = raw_input("* NEW FULL NAME = %s (OK?)" % full_name) # MODIFY IN FILE NAMES AND BARCODES JSON
				#if ui == '':
				print "* NEW FULL NAME = %s" % full_name
				barcodes_json[barcode]['sample'] = full_name
				# CHANGE SAMPLE FOLDER NAME
				subprocess.call(['mv',og_sample_name,full_name],cwd=run_folder)
				# CHANGE SAMPLE FILE NAMES
				os.chdir('%s/%s' % (run_folder,full_name))
				fl = glob.glob('*%s*.*' % (og_sample_name))
				for f in fl:
					newf = f.replace(og_sample_name,full_name)
					subprocess.call(['mv',f,newf],cwd='%s/%s' % (run_folder,full_name))
				modifications = True
				# SCRIPT ALAMUT
				vbs_scripts = glob.glob('*.vbs')
				for script in vbs_scripts:
					with open(script,'r') as s:
						lines = []
						for line in s:
							line = line.replace(og_sample_name,full_name)
							lines.append(line)
					with open(script,'w') as s:
						for line in lines:
							s.write(line)

	if modifications == True:					
		print '- WRITING BARCODES JSON...' # ouvrir le fichier en lecture open(...,'a') pour lire et modifier? pour eviter fichier temp et copie
		with open(run_folder+'/barcodes.json', 'w') as g:
			json.dump(barcodes_json,g,indent=4)		
