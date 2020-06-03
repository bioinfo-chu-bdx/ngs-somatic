#!/usr/bin/bash

echo "- git pull"
git pull origin master
echo "- download big files from gdrive"
echo "    - VariantBase.db"
$NGS_PIPELINE_BX_DIR/utility/gdrive-linux-x64 --service-account gdrive-key.json download '16VGtVtq5weh2AFfj-SXmpBBvBUiV1xwU' --path $NGS_PIPELINE_BX_DIR/variantBase/ --force
echo "    - cosmic.db"
$NGS_PIPELINE_BX_DIR/utility/gdrive-linux-x64 --service-account gdrive-key.json download '1v3-gIBN-1Us-LU8g5jtlu8huMZGc7Ybs' --path $NGS_PIPELINE_BX_DIR/variantBase/ --force
echo "    - variantList_ALL.json"
$NGS_PIPELINE_BX_DIR/utility/gdrive-linux-x64 --service-account gdrive-key.json download '1KecvpkoeMzSGM7gtrGmVwvzVEdhddMj_' --path $NGS_PIPELINE_BX_DIR/variantBase/variantList/ --force

# TODO !
# check if essentials programs and larges files (like databases) are installed (and pull / download from gdrive?)
