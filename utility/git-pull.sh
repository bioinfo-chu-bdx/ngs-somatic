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

echo "- git pull"
git pull origin master
echo "- download big files from gdrive"
$NGS_PIPELINE_BX_DIR/gdrive-linux-x64 --service-account gdrive-key.json download $NGS_PIPELINE_BX_DIR/variantBase/VariantBase.db
$NGS_PIPELINE_BX_DIR/gdrive-linux-x64 --service-account gdrive-key.json download $NGS_PIPELINE_BX_DIR/variantBase/cosmic.db
$NGS_PIPELINE_BX_DIR/gdrive-linux-x64 --service-account gdrive-key.json download $NGS_PIPELINE_BX_DIR/variantBase/variantList/variantList_ALL.json

# TODO !
# 2 : check if essentials programs and larges files (like databases) are installed (and pull / download from gdrive?)
