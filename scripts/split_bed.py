#!/usr/bin/env python
import os
import sys
import csv
from optparse import OptionParser

# This script split a bed for multi-threading purpose
# It should make clusters with approx. same size (in total base)
# BUT It prioritize keeping together (in same cluster) regions overlapping or close

### GATHERING PARAMETERS ############################################################

FNULL = open(os.devnull, 'w')
parser = OptionParser()
parser.add_option('-b', '--bed',			help="bed file",dest='bed') 
parser.add_option('-s', '--scatter-count',	help="scatter count",dest='scatter_count') 
parser.add_option('-o', '--output-folder',	help="output folder", dest='output_folder')
parser.add_option('-a', '--adjust',			help="adjust", dest='adjust',default=0)
parser.add_option('-m', '--chunk-spacing',	help="chunk-spacing", dest='chunk_spacing',default=250)
(options, args) = parser.parse_args()

min_between_chunk = int(options.chunk_spacing)
bed = options.bed
bed_name = bed.split('/')[-1].replace('.bed','')
scatter_count = int(options.scatter_count)
adjust = int(options.adjust)
bed_file = open(bed,'r')
bed_reader = csv.reader(bed_file,delimiter='\t')

chunks = {1:[]}

bed_total_base = 0
previous_end = 0
previous_chrm = 0
# COMPUTE BED TOTAL BASES
for line in bed_reader :
	if line[0].startswith('track'):
		continue
	chrm = line[0]
	start = int(line[1])
	end = int(line[2])
	interval_size = end - start
	if chrm == previous_chrm :
		if start <= previous_end :
			interval_size = end - previous_end
	previous_end = end
	previous_chrm = chrm
	bed_total_base += interval_size

print "\nbed_total_base = %s" % bed_total_base
chunk_size = (bed_total_base/scatter_count) + adjust
print "mean_chunk_size = %s" % chunk_size
# print "\nchunks (seek 12 for 12 threads):"
print "- splitting..."

bed_file.seek(0)
chunk_total_base = 0
previous_end = 0
previous_chrm = 0
chunk = 1
for line in bed_reader :
	if line[0].startswith('track'):
		continue
	chrm = line[0]
	start = int(line[1])
	end = int(line[2])
	if chrm == previous_chrm:
		if start <= previous_end : # SI chevauchement
			interval_size = end - previous_end
			chunk_total_base += interval_size
			# chunks[chunk].append([chrm,start,end])
			chunks[chunk].append(line)
			previous_end = end
			continue
		if (start-previous_end) < min_between_chunk : # SI espacement insuffisant entre interval
			interval_size = end - start
			chunk_total_base += interval_size
			# chunks[chunk].append([chrm,start,end])
			chunks[chunk].append(line)
			previous_end = end
			continue
	interval_size = end - start
	if ((chunk_total_base > chunk_size) or ((chunk_total_base+interval_size) > chunk_size)) and chunk<scatter_count:
		print "\t-chunk %s : size = %s" % (chunk,chunk_total_base)
		# next chunk
		chunk += 1
		chunks[chunk] = []
		# chunks[chunk].append([chrm,start,end])
		chunks[chunk].append(line)
		chunk_total_base = interval_size
	else:
		chunk_total_base += interval_size
		# chunks[chunk].append([chrm,start,end])
		chunks[chunk].append(line)
	previous_chrm = chrm
	previous_end = end

print "\t-chunk %s : size = %s" % (chunk,chunk_total_base)
i = 0
for chunk in chunks :
	chunk_file_name = '%04d-scattered.bed' % i
	chunk_file = open('%s/%s' % (options.output_folder,chunk_file_name),'w')
	chunk_writer = csv.writer(chunk_file,delimiter='\t')
	for line in chunks[chunk]:
		chunk_writer.writerow(line)
	chunk_file.close()
	i+=1
print "\n- chunks files writted in %s" % options.output_folder
