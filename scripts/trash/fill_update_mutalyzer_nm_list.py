#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import json
import csv
import urllib2
import xlrd
import openpyxl
import copy
import glob
import ast
import time
from subprocess import *
from optparse import OptionParser
    
parser = OptionParser()
parser.add_option('-b', '--bed',		help="Add new NM from annotated bed", dest='bed')
parser.add_option('-u', '--update',		help="Update mutalyzer_nm_list json : search new versions if existant", action='store_true' , default=False ,dest='update')
(options, args) = parser.parse_args()

jsn = open('/DATA/work/finalReport/mutalyzer_nm_list.json','r')
mutalyzer_nm_list = json.load(jsn)
jsn.close()

if not options.bed and not options.update:
	print "no options specified (--bed or/and --update)"
	
if options.bed:
	bed_file = open(options.bed,'r')
	bed_lines = bed_file.readlines()
	for i in range(1,len(bed_lines)):
		bed_line = bed_lines[i].split('\t')
		nm = bed_line[7].split('SUBMITTED_REGION=')[-1].split(';')[0].replace('\n','')
		if not (nm in mutalyzer_nm_list.keys()):
			print "* new nm : " + nm
			cp = 'c.20C>G'
			url_mutalyzer = "https://mutalyzer.nl/json/runMutalyzerLight?variant=%s:%s" % (nm,cp)
			f = urllib2.urlopen(url_mutalyzer)
			d = json.loads(f.read())
			version = d['sourceVersion']
			mutalyzer_nm_list[nm] = version
				
if options.update : 
	for nm in mutalyzer_nm_list:
		cp = 'c.20C>G'
		url_mutalyzer = "https://mutalyzer.nl/json/runMutalyzerLight?variant=%s:%s" % (nm,cp)
		f = urllib2.urlopen(url_mutalyzer)
		d = json.loads(f.read())
		version = d['sourceVersion']
		if mutalyzer_nm_list[nm] != version:
			print "* new version for " + nm
			print "\t %s -> %s" % (mutalyzer_nm_list[nm],version)
			mutalyzer_nm_list[nm] = version

jsn = open('/DATA/work/finalReport/mutalyzer_nm_list.json','w')		
json_text = json.dumps(mutalyzer_nm_list,indent=4)
jsn.write(json_text)
jsn.close()
