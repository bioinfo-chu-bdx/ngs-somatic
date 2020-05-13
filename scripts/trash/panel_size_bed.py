#!/usr/bin/env python
import sys

bedpath = sys.argv[1]
bed = open(bedpath,'r')
bed.next()

total_size = 0

previous_start = False
previous_stop = False
for line in bed:
	line = line.split('\t')
	start = int(line[1])
	stop = int(line[2])
	
	size = stop - start
	
	if previous_start:
		if previous_stop > start > previous_start :
			overlap = previous_stop - start
			size = size - overlap
		
	previous_start = start
	previous_stop = stop
	
	total_size = total_size + size

print "-" + bedpath.split('.bed')[0].split('/')[-1] + " : total size (without overlap) = " + str(total_size/1000.0) + " kb"
