##!/usr/bin/env python
import subprocess
import json
import sys

def json_loads_byteified(json_text):
    return _byteify(
        json.loads(json_text, object_hook=_byteify),
        ignore_dicts=True
    )

def _byteify(data, ignore_dicts = False):
    # if this is a unicode string, return its string representation
    if isinstance(data, unicode):
        return data.encode('utf-8')
    # if this is a list of values, return list of byteified values
    if isinstance(data, list):
        return [ _byteify(item, ignore_dicts=True) for item in data ]
    # if this is a dictionary, return dictionary of byteified keys and values
    # but only if we haven't already byteified it
    if isinstance(data, dict) and not ignore_dicts:
        return {
            _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True)
            for key, value in data.iteritems()
        }
    # if it's anything else, return it in its original form
    return data

results_file = open('/media/n06lbth/sauvegardes_pgm/FusionPlex_CTL/all_run_statistics.tsv','w')
results_file.write('RUN\tPATIENT\tTOTAL FRAGMENTS\tMAPPED READS\tDNA READS\tPERCENT\tRNA READS\tPERCENT\tRNA MEDIAN FRAGMENT LENGTH\n')

r = subprocess.check_output('curl -X GET "http://10.67.63.93/rest_api/jobs/" -H "accept: application/json" -H "authorization: Basic c2J0QGNodS1ib3JkZWF1eC5mcjpzYnQxMjNTQlQh" -H "X-CSRFToken: pBJkg5aAfS0LAb33LLQBwbVlnyBfIM0yzKsCLxw2xdVbzQaEGcFXNR31qdO00SIA"',shell=True)
j = json_loads_byteified(r)
for result in j['results']:
	job_name = result['name']
	job_id = result['job_id']
	r2 = subprocess.check_output('curl -X GET "http://10.67.63.93/rest_api/jobs/%s/" -H "accept: application/json" -H "authorization: Basic c2J0QGNodS1ib3JkZWF1eC5mcjpzYnQxMjNTQlQh" -H "X-CSRFToken: pBJkg5aAfS0LAb33LLQBwbVlnyBfIM0yzKsCLxw2xdVbzQaEGcFXNR31qdO00SIA"' % job_id,shell=True)
	j2 = json_loads_byteified(r2)
	for sample in j2['samples']:
		sample_name = sample['name'].replace('.unaligned','')
		sample_id = sample['id']
		r3 = subprocess.check_output('curl -X GET "http://10.67.63.93/rest_api/jobs/%s/samples/%s/results/read-stats" -H "accept: application/json" -H "authorization: Basic c2J0QGNodS1ib3JkZWF1eC5mcjpzYnQxMjNTQlQh" -H "X-CSRFToken: pBJkg5aAfS0LAb33LLQBwbVlnyBfIM0yzKsCLxw2xdVbzQaEGcFXNR31qdO00SIA"' % (job_id,sample_id),shell=True)
		j3 = json_loads_byteified(r3)
		total_fragments = j3['molbar_stats']['total_molbar_reads']
		for rstat in j3['read_stat_types']:
			if rstat['read_type'] == 'All Fragments':
				mapped_reads = rstat['mapped_num']
				break
		rna_median_fragment_length = j3['fragment_means_and_medians']['rna_median_length']
		for stat in j3['total_stats']:
			if stat['frag_type'] == 'All Fragments':
				adn_percent = stat['dna_reads_percent']
				adn_reads = stat['dna_reads']
				arn_percent = stat['rna_reads_percent']
				arn_reads = stat['rna_reads']
				break
		results_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (job_name,sample_name,total_fragments,mapped_reads,adn_reads,adn_percent,arn_reads,arn_percent,rna_median_fragment_length))

results_file.close()
# print json.dumps(j,indent=4)
