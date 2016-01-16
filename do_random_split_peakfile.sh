#!/bin/bash
#this converts 

#script adapted from here
#https://sites.google.com/site/anshulkundaje/projects/idr#TOC-IDR-ANALYSIS-ON-SELF-PSEUDOREPLICATES

#Randomly split the mapped reads in a tagAlign file into 2 parts (pseudoReplicates) with approximately equal number of reads
#INPUT gzipped tagAlign file
#OUTPUT: two random partitions of the original tagAlign, containing roughly the same number of reads
if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <FILENAME> <PATH>"
        echo "FILENAME - input tagAlign.gz filename"
	echo "PATH - where to store the output file - NO trailing backslash"
        exit
fi

filename=$1;
FILE_ID=`basename $1 ".tagAlign.gz"`;
#outputDir='/mapped/selfPseudoReps' # output directory for pseudoreplicate files
#outputStub='chipSampleRep1' # prefix name for pseudoReplicate files

nlines=$( zcat ${filename} | wc -l ) # Number of reads in the tagAlign file
nlines=$(( (nlines + 1) / 2 )); # half that number
zcat "${filename}" | shuf | split -d -l ${nlines} - "$2/${FILE_ID}" # This will shuffle the lines in the file and split it into two parts
gzip "$2/${FILE_ID}00"
gzip "$2/${FILE_ID}01"
mv "$2/${FILE_ID}00.gz" "$2/${FILE_ID}.pr1.tagAlign.gz"
mv "$2/${FILE_ID}01.gz" "$2/${FILE_ID}.pr2.tagAlign.gz"












#fileName='chipSampleRep1.tagAlign.gz' # input tagAlign file name
#outputDir='/mapped/selfPseudoReps' # output directory for pseudoreplicate files
#outputStub='chipSampleRep1' # prefix name for pseudoReplicate files

#nlines=$( zcat ${fileName} | wc -l ) # Number of reads in the tagAlign file
#nlines=$(( (nlines + 1) / 2 )) # half that number
#zcat "${fileName}" | shuf | split -d -l ${nlines} - "${outputDir}/${outputStub}" # This will shuffle the lines in the file and split it into two parts
#gzip "${outputDir}/${outputStub}00"
#gzip "${outputDir}/${outputStub}01"
#mv "${outputDir}/${outputStub}00.gz" "${outputDir}/${outputStub}.pr1.tagAlign.gz"
#mv "${outputDir}/${outputStub}01.gz" "${outputDir}/${outputStub}.pr2.tagAlign.gz"

