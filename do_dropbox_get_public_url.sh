#!/bin/bash

#2/6/2014
#dropbox command line
#get public url for all files in the dropbox public directory having a certain extension
#useful for UCSC

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <EXT>"
        echo "EXT - extension to consider when producing the URLs"
	echo "THIS ASSUMES dropbox is RUNNING (python ~/local/bin/dropbox.py start)"
	echo "THIS ASSUMES the dropbox public is /home/giuseppe/Dropbox/Public"
        exit
fi

PDATA='/home/giuseppe/Dropbox/Public';
PCODE='/home/giuseppe/local/bin'
EXT=$1

for FILE in ${PDATA}/*.${EXT};
	do
	source activate
	python ${PCODE}/dropbox.py puburl ${FILE}
done
