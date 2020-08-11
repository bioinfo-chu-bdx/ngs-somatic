#!/usr/bin/bash
run_path=$1
/usr/local/bin/multiqc -f -c $NGS_PIPELINE_BX_DIR/multiqc/custom_config.yaml -o $run_path/_multiqc -m 'fastqc' -m 'samtools' -m 'mosdepth' $run_path
#-m 'qualimap' -m 'samtools'
