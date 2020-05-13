#!/usr/bin/env python
import os
import sys
import json
import csv
import subprocess
import glob
import numpy
import matplotlib
matplotlib.use('Agg')
import pysam
from pylab import *
from optparse import OptionParser
import urllib2

sbt_sensitivity = open('/media/stuff/database_hotspot_sensitivity_27-09-2018.updated.tsv','r')
sbt_sensi_reader = csv.reader(sbt_sensitivity,delimiter='\t')

sbt_sensitivity_updated = open('/media/stuff/database_hotspot_sensitivity_01-10-2018.updated.tsv','w')
sbt_sensitivity_writer = csv.writer(sbt_sensitivity_updated,delimiter='\t')

with open('/DATA/work/global_parameters.json', 'r') as g:
	global_param = json.load(g)
with open(global_param['NM_data'], 'r') as nmdata:
	NM_version = json.load(nmdata)

sbt_sensi_reader.next()
sbt_sensitivity_writer.writerow(['Chr','Pos','Ref/Alt','ID','NM','Gene','Exon','c. (HGVS)','p. (HGVS)','Sensibilite'])
for line in sbt_sensi_reader:
	nm = line[4] + '.' + NM_version[line[4]]['version']
	cpos = line[7]
	cpos = cpos.replace('+','%2B')  # e.g : c.1959+1G>A
	cpos = cpos.replace('-','%2D')  # e.g : c.1912-2A>C
	url_mutalyzer = "https://mutalyzer.nl/json/runMutalyzer?variant=%s:%s" % (nm,cpos)
	f = urllib2.urlopen(url_mutalyzer)
	d = json.loads(f.read())

	try:
		visu = d['rawVariants'][0]['visualisation']
		
		visu1 = visu.split('\n')[0]
		leftseq1 = visu1.split(' ')[0][-3:]
		middleseq1 = visu1.split(' ')[1]
		rightseq1 = visu1.split(' ')[2][:3]
		
		visu2 = visu.split('\n')[1]
		leftseq2 = visu2.split(' ')[0][-3:]
		middleseq2 = visu2.split(' ')[1]
		rightseq2 = visu2.split(' ')[2][:3]
	
		print "%s:%s \t %s/%s" % (nm,cpos,middleseq1,middleseq2)
		sbt_sensitivity_writer.writerow([line[0],line[1],middleseq1+'/'+middleseq2,line[3],line[4],line[5],line[6],line[7],line[8],line[9]])
	except:
		errormessage = d['messages'][0]['message']
		print "%s:%s \t %s" % (nm,cpos,errormessage)
		sbt_sensitivity_writer.writerow([line[0],line[1],errormessage,line[3],line[4],line[5],line[6],line[7],line[8],line[9]])
