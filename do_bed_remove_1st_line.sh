#!/bin/bash

#mass remove 1st line from all the .bed files or other files in the same directory
#needed to clean up bed files produced by gem/gps

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <EXTENSION>"
        echo "You must specify 1) data path (e.g. /home/me/files) and 2) file extension (e.g. bed)"
        exit
fi

PDATA=$1;
EXT=$2;

for FILE in ${PDATA}/*.${EXT};
	do ID=`basename ${FILE} ".${EXT}"`;
	awk '{if (NR!=1) {print}}' ${FILE} > ${PDATA}/${ID}_pp.bed
done
