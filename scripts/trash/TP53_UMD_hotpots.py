#!/usr/bin/env python
import pysam

def reverse(seq):
	seq = seq.replace('A','W').replace('T','X').replace('G','Y').replace('C','Z')
	seq = seq.replace('W','T').replace('X','A').replace('Y','C').replace('Z','G')
	seq = seq[::-1]
	return seq

#hotspot_path = '/DATA/work/reference_files/Hotspots_TP53.bed'
UMD_path = '/media/stuff/UMD_variants_EU_Pathogenix_filtered.txt'

#hotspots = []
#with open(hotspot_path,'r') as hotbed:
	#hotbed.next()
	#for line in hotbed:
		#l = line.split('\t')
		#h = (l[0],l[1],l[2],l[6].split(';ANCHOR')[0])
		#hotspots.append(h)
		
with open(UMD_path,'r') as UMD:
	UMD.next()
	for line in UMD:
		l = line.split('\t')
		ref = reverse(l[5])
		obs = reverse(l[6])
		hid = l[2]
		if hid == '':
			hid = l[1]
		anchor = pysam.faidx('/DATA/work/hg19.fasta','%s:%s-%s' % ('chr17',str(int(l[3])-1),str(int(l[3])-1))).split('\n')[1]
		umd_hot = 'chr17\t%s\t%s\t%s\t0\t+\tREF=%s;OBS=%s;ANCHOR=%s\tTP53' % (str(int(l[3])-1),l[3],hid,ref,obs,anchor)
		
		print umd_hot
		#if umd_hot not in hotspots:
			#print umd_hot
