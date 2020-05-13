#!/usr/bin/env python
from operator import itemgetter 
import sqlite3
import json

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)
db_path = global_param['VariantBase']
db_con = sqlite3.connect(db_path)
db_con.row_factory = dict_factory
db_cur = db_con.cursor()


#db_cur.execute("SELECT geneName FROM Gene")
#db_gene = db_cur.fetchone()
#print db_gene['geneName']

#db_cur.execute("SELECT geneName FROM Gene")
#db_genes = db_cur.fetchall()
#for db_gene in db_genes:
	#print db_gene['geneName']


#variants=[
	#{
		#'chromosome':7,
		#'start':1242,
		#'stop':1248
	#},
	#{
		#'chromosome':1,
		#'start':184661,
		#'stop':1645674
	#},
	#{
		#'chromosome':11,
		#'start':4274272,
		#'stop':16455425227272674
	#},
	#{
		#'chromosome':7,
		#'start':2222,
		#'stop':2223
	#}
#]

#sorted_variants = sorted(sorted(sorted(variants,key=itemgetter('stop')),key=itemgetter('start')),key=itemgetter('chromosome'))

##sorted_variants = sorted(sorted(sorted(variants.items(), key=lambda v : int(v[2])), key=lambda v : int(v[1])), key=lambda v : int(v[0].replace('chr','').replace('X','23').replace('Y','24')))
#print sorted_variants

db_cur.execute("SELECT * FROM Variant WHERE variantID='chr4:1807922-1807922:G>A'")
db_variant = db_cur.fetchone()
print type(db_variant['gnomad'])
print db_variant['gnomad']

