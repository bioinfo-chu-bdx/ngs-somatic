#!/bin/bash
echo "- replace actual base with empty one"
cp /$NGS_PIPELINE_BX_DIR/variantBase/VariantBase.empty.db /$NGS_PIPELINE_BX_DIR/variantBase/VariantBase.db
echo "- insert_all_panel_gene_variantbase.py"
time python /$NGS_PIPELINE_BX_DIR/variantBase/insert_all_panel_gene_variantbase.py > /$NGS_PIPELINE_BX_DIR/variantBase/reannotation/logs/insert_all_panel_gene_variantbase.log.txt
echo "- insert_all_run_all_samples_variantbase.py"
time python /$NGS_PIPELINE_BX_DIR/variantBase/insert_all_run_all_samples_variantbase.py > /$NGS_PIPELINE_BX_DIR/variantBase/reannotation/logs/insert_all_run_all_samples_variantbase.log.txt
echo "- annotate_variantbase.step0.py --full"
time python /$NGS_PIPELINE_BX_DIR/variantAnnotation/annotate_variantbase.step0.py --full > /$NGS_PIPELINE_BX_DIR/variantBase/reannotation/logs/reannotate_all_variantbase.step0.log.txt
echo "- annotate_variantbase.step1.py --full"
time python /$NGS_PIPELINE_BX_DIR/variantAnnotation/annotate_variantbase.step1.py --full --output-folder /$NGS_PIPELINE_BX_DIR/variantBase/reannotation > /$NGS_PIPELINE_BX_DIR/variantBase/reannotation/logs/reannotate_all_variantbase.step1.log.txt
echo "- annotate_variantbase.step2.py"
time python /$NGS_PIPELINE_BX_DIR/variantAnnotation/annotate_variantbase.step2.py --annovar-results /$NGS_PIPELINE_BX_DIR/variantBase/reannotation/annovar/annovar.hg19_multianno.txt --vep-results /$NGS_PIPELINE_BX_DIR/variantBase/reannotation/vep/vep_output.tsv > /$NGS_PIPELINE_BX_DIR/variantBase/reannotation/logs/reannotate_all_variantbase.step2.log.txt
echo "- find_sample_pathology_SBT.py"
time python /$NGS_PIPELINE_BX_DIR/variantBase/find_sample_pathology_SBT.py > /$NGS_PIPELINE_BX_DIR/variantBase/reannotation/logs/find_sample_pathology_SBT.log.txt
echo "done."
