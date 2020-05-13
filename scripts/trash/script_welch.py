#!/usr/bin/env python
from scipy import stats
import csv

protfile_path = "/media/stuff/database_proteo.csv"
protfile = open(protfile_path,'r')
protfile_reader = csv.reader(protfile, delimiter = ';')

protfile_reader.next()
for line in protfile_reader:
	mr_list = []
	br_list = []
	protname = line[0]
	for i in range(2,14):
		mr_list.append(float(line[i]))
	for i in range(14,26):
		br_list.append(float(line[i]))
	
	statistic, pvalue = stats.ttest_ind(mr_list,br_list,equal_var=False)
	print protname+'\t'+str(statistic)+'\t'+str(pvalue)
