#!/usr/bin/python
import json
import os
import csv
import sqlite3

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
with open('%s/global_parameters.json' % pipeline_folder, 'r') as g:
	global_param = json.load(g)

db_path = global_param['cosmicDB']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()

cosmic_file = open('%s/variantAnnotation/annovar/humandb/hg19_cosmic90.txt' % pipeline_folder,'r')
cosmic_reader = csv.reader(cosmic_file,delimiter='\t')

for line in cosmic_reader:
	if line[0] == 'MT':
		continue
	chrm = 'chr%s' % line[0]
	start = int(line[1])
	stop = int(line[2])
	ref = line[3]
	alt = line[4]
	cosmicID = line[5].split(';')[0].split('ID=')[-1]
	db_cur.execute("INSERT INTO Cosmic (cosmicID, chr, start, stop, ref, alt) VALUES ('%s', '%s', %s, %s, '%s', '%s')" % (cosmicID,chrm,start,stop,ref,alt))
	
	tissues_occurences = line[5].split('OCCURENCE=')[-1]
	tissues = tissues_occurences.split(',')
	for t in tissues:
		tissue = t.split('(')[-1].split(')')[0]
		occurence = int(t.split('(')[0])
		db_cur.execute("UPDATE Cosmic SET '%s'=%s WHERE cosmicID='%s'" % (tissue,occurence,cosmicID))

db_con.commit()
db_con.close()
