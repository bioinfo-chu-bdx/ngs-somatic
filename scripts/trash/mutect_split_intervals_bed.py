#!/usr/bin/env python
import sys
import csv

bed = sys.argv[1]
num_chunks = int(sys.argv[2])
adjust = int(sys.argv[3])
min_between_chunk = 250

bed_file = open(bed,'r')
bed_reader = csv.reader(bed_file,delimiter='\t')

chunks = {1:[]}

# first compute total length of the bed
bed_reader.next()
bed_total_base = 0
previous_end = 0
previous_chrm = 0
for line in bed_reader :
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

print "bed_total_base = %s" % bed_total_base
chunk_size = (bed_total_base / num_chunks) - adjust
print "mean_chunk_size = %s" % chunk_size

bed_file.seek(0)
bed_reader.next()
chunk_total_base = 0
previous_end = 0
previous_chrm = 0
chunk = 1
for line in bed_reader :
	chrm = line[0]
	start = int(line[1])
	end = int(line[2])

	if chrm == previous_chrm:
		if start <= previous_end : # SI chevauchement
			interval_size = end - previous_end
			chunk_total_base += interval_size
			chunks[chunk].append([chrm,start,end])
			previous_end = end
			continue
			
		if (start-previous_end) < min_between_chunk : # SI espacement insuffisant entre interval
			interval_size = end - start
			chunk_total_base += interval_size
			chunks[chunk].append([chrm,start,end])
			previous_end = end
			continue
	
	interval_size = end - start
			
	if (chunk_total_base > chunk_size) or ((chunk_total_base+interval_size) > chunk_size):
		print "chunk %s : size = %s" % (chunk,chunk_total_base)
		# next chunk
		chunk += 1
		chunks[chunk] = []
		chunks[chunk].append([chrm,start,end])
		chunk_total_base = interval_size
	else:
		chunk_total_base += interval_size
		chunks[chunk].append([chrm,start,end])
		
	previous_chrm = chrm
	previous_end = end
		
print "chunk %s : size = %s" % (chunk,chunk_total_base)
	
l = 0		
for chunk in chunks :
	chunk_file_name = bed.split('/')[-1].replace('.bed','.chunk%s.intervals'%chunk)
	chunk_file = open('/DATA/work/reference_files/mutect/%s'%chunk_file_name,'w')
	chunk_writer = csv.writer(chunk_file)
	for interval in chunks[chunk]:
		itv = "%s:%s-%s" % (interval[0],interval[1],interval[2])
		chunk_writer.writerow([itv])
	chunk_file.close()
	
