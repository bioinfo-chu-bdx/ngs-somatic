date
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -geneanno -buildver hg19 -dbtype refGene -outfile /DATA/work/scripts/tests/variants.annotated.refGene -exonsort /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ --hgvs &
echo "1"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype cosmic86 -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ &
echo "2"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype avsnp150 -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ &
echo "3"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype intervar_20180118 -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ -otherinfo &
echo "4"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype clinvar_20180603 -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ -otherinfo &
echo "5"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype nci60 -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ &
echo "6"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype esp6500siv2_all -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ &
echo "7"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ &
echo "8"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_eur -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ &
echo "9"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_amr -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ &
echo "10"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_afr -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ &
echo "11"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_eas -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ &
echo "12"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_sas -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ &
echo "13"
perl /DATA/work/variantAnnotation/annovar/annotate_variation.pl -filter -dbtype dbnsfp33a -buildver hg19 -outfile /DATA/work/scripts/tests/variants.annotated /DATA/work/scripts/tests/annovar_input.tsv /DATA/work/variantAnnotation/annovar/humandb/ -otherinfo &
wait
date
