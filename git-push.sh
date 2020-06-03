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

echo "git commit ...";
git commit -m "$USER-`date '+%Y-%m-%d'`";
echo "git push ...";
git push -f origin master;
expect "Username for 'https://github.com': "
send "thomas.bandres@chu-bordeaux.fr"
expect "Password for 'https://thomas.bandres@chu-bordeaux.fr@github.com': "
send "8bc%VXMOtZA9"

# to pull : git pull origin master
