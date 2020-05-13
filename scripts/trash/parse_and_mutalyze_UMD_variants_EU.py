#!/usr/bin/env python
import os
import json
import urllib2
import csv

def request_mutalyzer(nm,cp):
    nuc_hgvs = ''
    prot_hgvs = ''
    failstring = ''
    prot_index = 0
    try:
        cp = cp.replace('+','%2B')  # e.g : c.1959+1G>A
        cp = cp.replace('-','%2D')  # e.g : c.1912-2A>C
        succes = True
        url_mutalyzer = "https://mutalyzer.nl/json/runMutalyzerLight?variant=%s:%s" % (nm,cp)
        f = urllib2.urlopen(url_mutalyzer)
        d = json.loads(f.read())
        if d['transcriptDescriptions']:
            if len(d['transcriptDescriptions']) == 1:
                nuc_description = d['transcriptDescriptions'][0]
                nuc_hgvs = nuc_description.split(':')[-1]
            else: # with LRG, several entries are possible 
                for i in range(len(d['transcriptDescriptions'])):
                    nuc_description = d['transcriptDescriptions'][i]
                    lrg = nuc_description.split(':')[0]
                    if lrg == nm:
                        nuc_hgvs = nuc_description.split(':')[-1]
                        prot_index = i
                        break
        else:
            succes = False
            print "%s:%s : Cannot get p. from Mutalyzer." % (nm,cp)
        if d['proteinDescriptions']:
            prot_description = d['proteinDescriptions'][prot_index]
            prot_hgvs = prot_description.split(':')[-1]
        else:
            succes = False
            print "%s:%s : Cannot get p. from Mutalyzer." % (nm,cp)
        if not succes:
            failstring = "Mutalyzer : "
            for message in d['messages']:
                if message['errorcode'] != 'WNOVER':
                    failstring = failstring + message['message']
    except Exception as e: 
		print "Unexpected error from request_mutalyzer(%s,%s) : %s" % (nm, cp, e)
		print "Failed to get data from Mutalyzer Name Checker"
        
    return (nuc_hgvs,prot_hgvs,failstring)

UMD_variants_EU = '/media/stuff/UMD_variants_EU.txt'
umdfile = open(UMD_variants_EU, 'r') 
csvreader = csv.reader(umdfile, delimiter = '\t')

results_path = '/DATA/work/finalReport/TP53_UMD_variants_EU.tsv'
results = open(results_path,'w')
csvwriter = csv.writer(results, delimiter = '\t')

csvreader.next()
nm = 'NM_000546.5'

csvwriter.writerow(['cDNA','Protein','Pathogenicity','Final comment'])
for line in csvreader:
	cpos = line[8]
	prot = line[24]
	final_comment = line[-1]
	pathogenicity = line[-2]
	if ('del' in cpos) and ('ins' in cpos):
		cpos = cpos.split('del')[0] + 'delins' + cpos.split('ins')[-1]
	elif ('del' in cpos):
		cpos = cpos.split('del')[0] + 'del'
	elif ('ins(' in cpos):
		csvwriter.writerow([cpos,prot,pathogenicity,final_comment,'hgvs not checked'])
		continue
	elif ('+' in cpos) or ('-' in cpos):
		csvwriter.writerow([cpos,prot,pathogenicity,final_comment,'hgvs not checked'])
		continue
	hgvs = request_mutalyzer(nm,cpos)
	if hgvs[0] :
		cdna = hgvs[0]
		prot = hgvs[1]
		csvwriter.writerow([cdna,prot,pathogenicity,final_comment])
	else:
		print "Mutalyzer failed for %s. Keeping original descriptions." % cpos
		print "-> %s \n" % hgvs[2]
		cdna = cpos
		prot = prot
		csvwriter.writerow([cdna,prot,pathogenicity,final_comment,hgvs[2]])



