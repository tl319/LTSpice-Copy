#!/bin/bash
# Usage: remove all utility bills pdf file password 

for f in *.net
do
	echo "processing -$f"
	cat $f | ./output
	SUBSTR=$(echo $f | cut -d'.' -f 1)
	SUBSTR+=".txt"
	cat $SUBSTR 
	printf "\n\n\n" 
done
