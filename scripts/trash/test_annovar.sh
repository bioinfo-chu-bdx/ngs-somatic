#!/usr/bin/env bash
vcf='/media/stuff/runs/test_3samples_SBT/ARISTIDE-AU453F/intermediate_files/tvc_only_hotspot/TSVC_variants.vcf'

input='/media/stuff/runs/test_3samples_SBT/ARISTIDE-AU453F/intermediate_files/annovar/annovar_input.tsv'
output='/media/stuff/runs/test_3samples_SBT/ARISTIDE-AU453F/intermediate_files/annovar/variants.test.annotated'

#time perl /DATA/work/variantAnnotation/annovar/table_annovar.pl $input /DATA/work/variantAnnotation/annovar/humandb/ -buildver 'hg19' -out $output -protocol 'refGeneWithVer,cosmic89,avsnp150,intervar_20180118,clinvar_20190305,nci60,esp6500siv2_all,1000g2015aug_all,gnomad211_genome,exac03,dbnsfp33a' -operation 'g,f,f,f,f,f,f,f,f,f,f' -argument '--hgvs,,,,,,,,,,' -nastring '.' -polish -xref /DATA/work/variantAnnotation/annovar/example/gene_xref.txt
#time perl /DATA/work/variantAnnotation/annovar/table_annovar.pl $input /DATA/work/variantAnnotation/annovar/humandb/ -buildver 'hg19' -out $output -protocol 'refGeneWithVer' -operation 'g' -argument '' -nastring '.' -polish -xref /DATA/work/variantAnnotation/annovar/example/gene_xref.txt


vepinput='/media/stuff/runs/test_3samples_SBT/ARISTIDE-AU453F/intermediate_files/vep/vep_input.tsv'
vepoutput='/media/stuff/runs/test_3samples_SBT/ARISTIDE-AU453F/intermediate_files/vep/vep_output2.tsv'

time perl /DATA/work/variantAnnotation/VEP/ensembl-vep/vep --offline --dir_cache /DATA/work/variantAnnotation/VEP/cache/ --force_overwrite --refseq --numbers --hgvs --hgvsg --no_escape --variant_class --sift 'p' --polyphen 'p' --af_1kg --freq_pop '1KG_ALL' --freq_pop 'gnomAD' --af_esp --symbol --tab --pubmed --fasta /DATA/work/variantAnnotation/VEP/cache/homo_sapiens/96_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --exclude_predicted --no_stats --verbose -i $vcf -o $vepoutput
