#!/usr/bin/env python
import os
import sys
import csv

csv.field_size_limit(sys.maxsize)

######
# PROTOCOLE POUR GENERER DES HOTSPOTS :
# 0 - Download from COSMIC :
# - CosmicMutantExport.tsv
# - CosmicCodingMuts.vcf
# 1 - utiliser ce script : python generate_hotspots.py <CosmicTSV> <CosmicVCF> <bedTarget> ou ajout manuel des variables
# 2 - tvcutils : variantCaller/bin/tvcutils prepare_hotspots -v hotspot.vcf -o output.vcf -r hg19.fasta -u target_bed


count_pass = 10

pipeline_folder = os.environ['NGS_PIPELINE_BX_DIR']
CosmicMutantExport_tsv_path = '/media/stuff/COSMIC/v92/CosmicMutantExport.tsv'#sys.argv[1]
CosmicCodingMuts_vcf_path = '/media/stuff/COSMIC/v92/CosmicCodingMuts.vcf'#sys.argv[2]
target_bed_path = '%s/reference_files/LAM2018_IAD143291_236_Designed_with_NM.bed' % pipeline_folder #sys.argv[3]
output = '%s/reference_files/annotate_bed/LAM2018_IAD143291_236.hotspots.10.vcf' % pipeline_folder


####################################################################################################################################################
temoin_hotspots = []
primary_sites = []
histology_subtypes = []
site_subtype_excluded = []

#### PRIMARY SITE AND HISTOLOGY SUBTYPE (UNCOMMENT IF NEEDED)
## SBT
#temoin_hotspots = ['COSM584','COSM1332933','COSM577','COSM564','COSM24850','COSM28056','COSM28055','COSM759','COSM760','COSM763','COSM1420865','COSM328026','COSM775','COSM12464','COSM715','COSM29446','COSM723','COSM716','COSM24842','COSM724','COSM721','COSM1428724','COSM719','COSM12418','COSM1430085','COSM22415','COSM1430086','COSM587613','COSM736','COSM28052','COSM24637','COSM41602','COSM1326','COSM96867','COSM96885','COSM1430136','COSM1275','COSM1290','COSM1299','COSM1304','COSM12706','COSM20402','COSM19194','COSM133767','COSM13177','COSM41905','COSM6239','COSM13979','COSM53194','COSM13182','COSM17570','COSM6223','COSM28603','COSM13190','COSM12986','COSM28610','COSM53291','COSM13424','COSM6227','COSM13430','COSM6224','COSM6213','COSM14070','COSM13008','COSM706','COSM710','COSM29633','COSM1447462','COSM697','COSM43064','COSM1214928','COSM1447471','COSM1330154','COSM700','COSM48565','COSM695','COSM691','COSM476','COSM471','COSM467','COSM462','COSM450','COSM27986','COSM1448625','COSM6262','COSM36912','COSM36903','COSM499','COSM495','COSM483','COSM19940','COSM554','COSM546','COSM14208','COSM521','COSM507','COSM1235478','COSM1235479','COSM20959','COSM21212','COSM21570','COSM27316','COSM20944','COSM25229','COSM29005','COSM21355','COSM21360']
#primary_sites = ['lung','large_intestine','skin']
#site_subtype_excluded = ['caecum','anorectal','anus','appendix']
## LAM
temoin_hotspots = ['COSM1411076','COSM1681610','COSM3732385','COSM1676499','COSM53042','COSM3259655','COSM783','COSM1418772','COSM28747','COSM33733','COSM24437','COSM12600','COSM532','COSM158604','COSM583','COSM1681955','COSM133120','COSM211643','COSM10812']
primary_sites = ['haematopoietic_and_lymphoid_tissue']
histology_subtypes = ['acute_leukaemia_of_ambiguous_lineage',
	'acute_leukaemic_transformation_of_primary_myelofibrosis',
	'acute_leukaemic_transformation_of_essential_thrombocythaemia',
	'acute_leukaemic_transformation_of_polycythaemia_vera',
	'acute_leukaemic_transformation_of_myeloproliferative_neoplasm',
	'acute_leukaemic_transformation_of_chronic_myelomonocytic_leukaemia',
	'acute_myeloid_leukaemia',
	'acute_myeloid_leukaemia_associated_with_MDS',
	'acute_myeloid_leukaemia_therapy_related',
	'acute_myeloid_leukaemia_myelodysplastic_syndrome_therapy_related_NOS',
	'myelodysplastic_syndrome',
	'myelodysplastic-myeloproliferative_neoplasm-unclassifiable',
	'myelodysplastic_syndrome_therapy_related',
	'myelodysplastic-myeloproliferative_neoplasm',
	'myeloproliferative_neoplasm',
	'myeloproliferative_neoplasm_unclassifiable'
	]
## ABL1	
#primary_sites = ['haematopoietic_and_lymphoid_tissue']
#histology_subtypes = ['chronic_myeloid_leukaemia',
	#'chronic_myelomonocytic_leukaemia',
	#'blast_phase_chronic_myeloid_leukaemia',
	#]
## FLT3		
#primary_sites = ['haematopoietic_and_lymphoid_tissue']
#histology_subtypes = ['acute_leukaemia_of_ambiguous_lineage',
	#'acute_leukaemic_transformation_of_primary_myelofibrosis',
	#'acute_leukaemic_transformation_of_essential_thrombocythaemia',
	#'acute_leukaemic_transformation_of_polycythaemia_vera',
	#'acute_leukaemic_transformation_of_myeloproliferative_neoplasm',
	#'acute_leukaemic_transformation_of_chronic_myelomonocytic_leukaemia',
	#'acute_myeloid_leukaemia',
	#'acute_myeloid_leukaemia_associated_with_MDS',
	#'acute_myeloid_leukaemia_therapy_related',
	#'acute_myeloid_leukaemia_myelodysplastic_syndrome_therapy_related_NOS',
	#'essential_thrombocythaemia',
	#'myelodysplastic_syndrome',
	#'myelodysplastic-myeloproliferative_neoplasm-unclassifiable',
	#'myelodysplastic_syndrome_therapy_related',
	#'myelodysplastic-myeloproliferative_neoplasm',
	#'myeloproliferative_neoplasm',
	#'myeloproliferative_neoplasm_unclassifiable'
	#]
####################################################################################################################################################

CosmicMutantExport_file = open(CosmicMutantExport_tsv_path,'r')
CosmicMutantExport_reader = csv.reader(CosmicMutantExport_file,delimiter='\t')

CosmicCodingMuts_file = open(CosmicCodingMuts_vcf_path,'r')
CosmicCodingMuts_reader = csv.reader(CosmicCodingMuts_file,delimiter='\t')

target_bed_file = open(target_bed_path,'r')
target_bed_reader = csv.reader(target_bed_file,delimiter='\t')

results_hotspot_vcf = open(output,'w')
results_hotspot_writer = csv.writer(results_hotspot_vcf, delimiter = '\t')

bed_pos = []
target_bed_reader.next()
for line in target_bed_reader:
	if line[0] == 'chrX':
		chrm = 23
	elif line[0] == 'chrY':
		chrm = 24
	else:
		chrm = int(line[0].replace('chr',''))
	start = int(line[1])
	end = int(line[2])
	bed_pos.append((chrm,start,end))

CosmicMutantExport_reader.next()
cosm2keep = {}
i=0
for line in CosmicMutantExport_reader:
	i=i+1
	gene_name = line[0]
	mutation_genomic_position = line[25]
	mutation_description = line[21]
	Mutation_CDS = line[19]
	primary_site = line[7]
	site_subtype = line[8]
	histology_subtype = line[12]
	mutation_id = line[16] # COSV
	legacy_mutation_id = line[17] # COSM and other shit
	mutation_somatic_status = line[29]
	
	if '_ENST' in gene_name:
		continue
	if mutation_genomic_position == '':
		continue
	chrm = int(mutation_genomic_position.split(':')[0])
	start = int(mutation_genomic_position.split(':')[-1].split('-')[0])
	end = int(mutation_genomic_position.split('-')[-1])
	pos_ok = False
	for pos in bed_pos:
		if (chrm == pos[0]) and ( (pos[2] > start > pos[1]) or (pos[2] > end > pos[1]) ):
			pos_ok = True
			break
	if not pos_ok:
		continue

	if legacy_mutation_id in temoin_hotspots:
		if mutation_id not in cosm2keep:
			cosm2keep[mutation_id] = 1
		else:
			cosm2keep[mutation_id] = cosm2keep[mutation_id] + 1
		continue

	if mutation_description in ['Substitution - coding silent','Whole gene deletion','Unknown','Deletion - In frame','Insertion - In frame','Complex - insertion inframe','Complex - deletion inframe']:
		continue
		
	if '?' in Mutation_CDS:
		continue

	if primary_site not in primary_sites:
		continue
		
	if site_subtype in site_subtype_excluded:
		continue
		
	if histology_subtype not in histology_subtypes:
		continue
		
	if mutation_somatic_status in ['Variant of unknown origin','Not specified']:
		continue
		
	if mutation_id not in cosm2keep:
		cosm2keep[mutation_id] = 1
	else:
		cosm2keep[mutation_id] = cosm2keep[mutation_id] + 1

print str(len(cosm2keep.keys())) + ' variants retained over ' + str(i) + ' lines'
final_cosm2keep = []

for cosm in cosm2keep:
	if cosm2keep[cosm] >= count_pass:
		final_cosm2keep.append(cosm)
print str(len(final_cosm2keep)) + ' variants retained with %s+ occurences' % count_pass

### WRITE OUTPUT
for line in CosmicCodingMuts_reader:
	if line[0].startswith('#'):
		continue
	if line[2] in final_cosm2keep:
		results_hotspot_writer.writerow(line)


CosmicMutantExport_file.close()
CosmicCodingMuts_file.close()
target_bed_file.close()
results_hotspot_vcf.close()
