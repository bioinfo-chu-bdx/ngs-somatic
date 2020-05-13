#!/usr/bin/env python
import json
import sys

param_json1 = sys.argv[1]
param_json2 = sys.argv[2]

with open(param_json1, 'r') as g:
	json1 = json.load(g)
with open(param_json2, 'r') as g:
	json2 = json.load(g)
	
print "Differences in VariantCallers configurations :"
print "Parameter Name\tParam1 : %s\tParam2 : %s" % (param_json1.split('/')[-1],param_json2.split('/')[-1])
for section in ['freebayes','torrent_variant_caller','long_indel_assembler']:
	print section
	for key in json1[section].keys():
		if key in json2[section].keys():
			if (str(json1[section][key]) != str(json2[section][key])):
				print "%s\t%s\t%s"% (key,str(json1[section][key]),str(json2[section][key]))
			del json1[section][key],json2[section][key]
		else:
			print "%s\t%s\t%s"% (key,str(json1[section][key]),'NA')
	for key in json2[section].keys():
		print "%s\t%s\t%s"% (key,'NA',str(json2[section][key]))
