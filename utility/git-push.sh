#!/usr/bin/bash

##########################################################
confirm() {
# syntax: confirm [<prompt>]
# Prompts the user to enter Yes or No and returns 0/1.

  local _prompt _default _response
 
  if [ "$1" ]; then _prompt="$1"; else _prompt="Do you wish to continue?"; fi
  _prompt="$_prompt [y/n] ?"
 
  # Loop forever until the user enters a valid response (Y/N or Yes/No).
  while true; do
    read -r -p "$_prompt " _response
    case "$_response" in
      [Yy][Ee][Ss]|[Yy]) # Yes or Y (case-insensitive).
        return 0
        ;;
      [Nn][Oo]|[Nn])  # No or N.
        exit
        ;;
      *) # Anything else (including a blank) is invalid.
        ;;
    esac
  done
}
##########################################################

echo "- git add ...";
git add -A;
echo "- git status : ";
git status;
confirm;
echo "- git commit ...";
git commit -m "$USER-`date '+%Y-%m-%d'`";
echo "- git push ...";
git push -f origin master;
# thomas.bandres@chu-bordeaux.fr
# 8bc%VXMOtZA9
echo "- upload big files to gdrive :"
echo "    - VariantBase.db"
$NGS_PIPELINE_BX_DIR/utility/gdrive-linux-x64 --service-account gdrive-key.json update '16VGtVtq5weh2AFfj-SXmpBBvBUiV1xwU' $NGS_PIPELINE_BX_DIR/variantBase/VariantBase.db
echo "    - cosmic.db"
$NGS_PIPELINE_BX_DIR/utility/gdrive-linux-x64 --service-account gdrive-key.json update '1v3-gIBN-1Us-LU8g5jtlu8huMZGc7Ybs' $NGS_PIPELINE_BX_DIR/variantBase/cosmic.db
echo "    - variantList_ALL.json"
$NGS_PIPELINE_BX_DIR/utility/gdrive-linux-x64 --service-account gdrive-key.json update '1KecvpkoeMzSGM7gtrGmVwvzVEdhddMj_' $NGS_PIPELINE_BX_DIR/variantBase/variantList/variantList_ALL.json
#echo "    - FalsePositives.csv" # fichier petit deja sur github pas besoin
#$NGS_PIPELINE_BX_DIR/utility/gdrive-linux-x64 --service-account gdrive-key.json update '1N98azA__fb3uE4AMkz2aczpfLc2LJcRY' $NGS_PIPELINE_BX_DIR/finalReport/FalsePositives.csv

echo "... done."

# Note : to avoir constant login/password ask :
# git config --global credential.helper cache
# git config --global credential.helper 'cache --timeout=100000000'
