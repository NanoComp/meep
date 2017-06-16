#!/bin/bash

###################################################
# FixWikiConversion.sh -- fix broken elements of
# markdown files produced by conversion from wiki
#
# Homer Reid 6/2017
###################################################

if [ $# -ne 1 ]
then
  echo "usage: FixWikiConversion.sh File.md"
  exit
fi

FILE=$1
/bin/cp ${FILE} /tmp/${FILE}

vim -e ${FILE} <<-EndOfHere
	:so Replacements2.vim
	:wq!
EndOfHere
