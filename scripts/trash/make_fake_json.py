#!/usr/bin/python
import sys
import os
import glob
import json

run_folder = sys.argv[1]
json_data = {}

dir_elements = os.listdir(run_folder)
for element in dir_elements:
	if not os.path.isdir(run_folder + '/' + element):
		continue
	bamfile = glob.glob('%s/%s/*.bam' % (run_folder,element))
	if bamfile:
		barcode = 'IonXpress_' + bamfile[0].split('IonXpress_')[-1].split('.bam')[0].split('.processed')[0]
		json_data[barcode] = {}
		json_data[barcode]['sample'] = element
		json_data[barcode]['target_region_filepath'] = '/fake/fake/unmerged/detail/IAD108862_231_Designed.with_NM.bed'#IAD78219_237_Designed_with_NM.bed'#IAD62716_182_Designed_with_NM.bed'#IAD37093_Designed_with_NM.bed'
		json_data[barcode]['hotspot_filepath'] = '/fake/fake/unmerged/detail/SBT_Hotspots_ColonLung_sept2015_v4.bed'
		json_data[barcode]['barcode_description'] = 'lib1'

jsonfile = open(run_folder + '/barcodes.json','w')
jsonfile.write(json.dumps(json_data, sort_keys=True,indent=4))
jsonfile.close
