#!/bin/bash

echo "- annotate_variantbase.step0.py --new"
time python $NGS_PIPELINE_BX_DIR/variantAnnotation/annotate_variantbase.step0.py --new | tee $NGS_PIPELINE_BX_DIR/variantBase/reannotation/logs/reannotate_all_variantbase.step0.log.txt
echo "- annotate_variantbase.step1.py --new"
time python $NGS_PIPELINE_BX_DIR/variantAnnotation/annotate_variantbase.step1.py --new --output-folder $NGS_PIPELINE_BX_DIR/variantBase/reannotation | tee /$NGS_PIPELINE_BX_DIR/variantBase/reannotation/logs/reannotate_all_variantbase.step1.log.txt
echo "- annotate_variantbase.step2.py"
time python $NGS_PIPELINE_BX_DIR/variantAnnotation/annotate_variantbase.step2.py --annovar-results $NGS_PIPELINE_BX_DIR/variantBase/reannotation/annovar/annovar.hg19_multianno.txt --vep-results $NGS_PIPELINE_BX_DIR/variantBase/reannotation/vep/vep_output.tsv | tee /$NGS_PIPELINE_BX_DIR/variantBase/reannotation/logs/reannotate_all_variantbase.step2.log.txt
echo "- annotate_variantbase.custom_annotation.py"
time python $NGS_PIPELINE_BX_DIR/variantAnnotation/annotate_variantbase.custom_annotation.py | tee $NGS_PIPELINE_BX_DIR/variantBase/reannotation/logs/reannotate_all_variantbase.only_custom_annotation.log.txt
echo "- find_sample_pathology_SBT.py"
time python $NGS_PIPELINE_BX_DIR/variantBase/find_sample_pathology_SBT.py | tee $NGS_PIPELINE_BX_DIR/variantBase/reannotation/logs/find_sample_pathology_SBT.log.txt
echo "done."
