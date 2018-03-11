#!/bin/bash
echo "Warning: files may be overwiten silently"
read -p "Continue (y/n)? " choice
case "$choice" in 
  y|Y ) for fn in data/PGP/*;
	do
	echo "unzipping sample ${fn}"
	unzip ${fn} -d data/PGP 
	done;;
  n|N ) echo "did nothing";;
  * ) echo "invalid";;
esac 
