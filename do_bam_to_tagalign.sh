#!/bin/bash

#from 
#https://sites.google.com/site/anshulkundaje/projects/idr#TOC-IDR-ANALYSIS-ON-SELF-PSEUDOREPLICATES
#INPUT: bam file from bwa, stampy or bowtie
#OUTPUT: tagAlign file for SPP or phantompeakcall or Anshul Kundaje's other stuff


if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <FILENAME.bam> <PATH>"
        echo "FILENAME.bam: file to process"
	echo "PATH: where to store the output file - NO trailing backslash"
        exit
fi

PCODE='/home/giuseppe/local/bin';
PBEDTOOLS='/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin'
FILE_ID=`basename $1 ".bam"`;

${PCODE}/samtools view -b -F 1548 -q 30 $1 | ${PBEDTOOLS}/bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > $2/${FILE_ID}.tagAlign.gz

