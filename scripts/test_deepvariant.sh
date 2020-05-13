sudo docker run \
  -v /home/t8illumy/ngs-pipeline-bx/reference_files/hg19/hg19.fasta \
  -v /home/t8illumy/ngs-pipeline-bx/reference_files/SureSelect-HEMATO-v5.sorted.annotated.bed \
  -v /media/stuff/runs/run_test_illumina/AcroMetrix/ \
  google/deepvariant:0.9.0 \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WES \
  --ref=hg19.fasta \
  --reads=AcroMetrix.bam \
  --regions=SureSelect-HEMATO-v5.sorted.annotated.bed \
  --output_vcf=intermediate_files/deepvariant/AcroMetrix.deepvariant.output.vcf.gz \
  --output_gvcf=intermediate_files/deepvariant/AcroMetrix.deepvariant.output.g.vcf.gz \
  --num_shards=12
