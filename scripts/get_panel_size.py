#!/usr/bin/python

#panel_path = '%s/reference_files/Target_ColonLung_v10_IAD172906_231.bed'
panel_path = '/media/stuff/fastq/SureSelect-HEMATO-v5.sorted.bed'
panel_file = open(panel_path,'r')

chr_count = {'chr1':[],'chr2':[],'chr3':[],'chr4':[],'chr5':[],'chr6':[],'chr7':[],'chr8':[],'chr9':[],'chr10':[],'chr11':[],'chr12':[],'chr13':[],'chr14':[],'chr15':[],'chr16':[],'chr17':[],'chr18':[],'chr19':[],'chr20':[],'chr21':[],'chr22':[],'chrX':[],'chrY':[]}
base_count = 0
for line in panel_file.readlines():
	if not line.startswith('chr'):
		continue
	line = line.split('\t')
	print line
	chromosome = line[0]
	start = int(line[1])
	stop = int(line[2])
	for position in range(start+1,stop+1):
		if position not in chr_count[chromosome]:
			chr_count[chromosome].append(position)
			base_count += 1
print "%s = %s pb" % (panel_path,base_count)

#import sqlite3
#import json

#def dict_factory(cursor, row):
    #d = {}
    #for idx, col in enumerate(cursor.description):
        #d[col[0]] = row[idx]
    #return d
    
    
#panel = 'Target_ColonLung_v10_IAD172906_231.bed'
    
	
#with open('%s/global_parameters.json', 'r') as g:
	#global_param = json.load(g)
#db_path = global_param['VariantBase']

#db_con = sqlite3.connect(db_path)
#db_con.row_factory = dict_factory
#db_cur = db_con.cursor()
	
#db_cur.execute("SELECT chromosome,start,stop FROM TargetedRegion INNER JOIN Panel ON TargetedRegion.panel = Panel.panelID WHERE panel='%s' ORDER BY start" % panel)
#db_target_regions = db_cur.fetchall()

#base_count = []
#for db_target_region in db_target_regions:
	#chromosome = db_target_region['chromosome']
	#start = db_target_region['start']
	#stop = db_target_region['stop']
	#for position in range(start+1,stop+1):
		#chrpos = '%s:%s' % (chromosome,position)
		#if chrpos not in base_count:
			#base_count.append(chrpos)

#print "%s = %s pb" % (panel,len(base_count))
