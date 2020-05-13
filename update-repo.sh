#!/bin/bash

#	-SOURCE = dossier a sauvegarder
#	-DEST = nom du dossier de destination

SOURCE=$NGS_PIPELINE_BX_DIR/
DEST=/home/t8illumy/repository/

rsync -rltsOv \
--exclude 'reference_files/hg19/*' \
--exclude 'reference_files/mutect/af-only-gnomad.raw.sites.b37.vcf.gz' \
--exclude 'reference_files/mutect/_old' \
--exclude 'variantAnnotation/annovar/humandb/*' \
--exclude 'variantAnnotation/VEP/cache/homo_sapiens/*' \
--exclude 'variantAnnotation/VEP/cache/homo_sapiens_refseq/*' \
--exclude 'variantAnnotation/VEP/ensembl-vep/t/testdata/*' \
--exclude 'variantAnnotation/VEP/ensembl-vep/t/test-genome-DBs/*' \
--exclude 'variantAnnotation/VEP/ensembl-vep/.git/*' \
--exclude 'variantAnnotation/temp/*' \
--exclude 'variantBase/VariantBase.db' \
--exclude 'variantBase/VariantBase.db.*.backup' \
--exclude 'variantBase/old_files/*' \
--exclude 'variantBase/reannotation/*' \
--exclude 'seqrepo/*' \
--exclude 'pindel/test' \
--exclude 'pindel/demo' \
--exclude 'gatk/gatk-package*' \
$SOURCE $DEST

# the exclude path is relative to the source path, NOT /
#--delete \

echo "- Updating local repos done"

tar -czvf '/home/t8illumy/repository/ngs-pipeline-bx.tar.gz' $DEST --exclude $DEST'.git' --exclude $DEST'ngs-pipeline-bx.tar.gz'

echo "- Generating archive tar.gz done"

# ensuite :
# cd /DATA/repository/ngs-pipeline-bx/
# git add -A
# git commit -m "message"
# git push -u origin master
