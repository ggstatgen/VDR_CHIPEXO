#!/bin/bash

#It outputs 2 numbers: average coverage for a coveraged base (sum/cnt) AND total base coverage (sum)
#/net/isi-scratch/giuseppe/tools/samtools/samtools depth /net/isi-scratch/giuseppe/VDR/MY_FASTQ_CUTADAPT/BAM_STAMPY/VDR_GM06986_Peconic20301_trimmed.sorted.bam | awk '{sum+=$3;cnt++}END{print sum/cnt" "sum}'

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "You must specify the data path for the bam files (e.g. /home/me/files)"
        exit
fi

PDATA=$1;
#location of fastqc executable, change if needed
PCODE='/net/isi-scratch/giuseppe/tools';


for FILE in ${PDATA}/*.bam;
	do
	${PCODE}/samtools/samtools depth ${FILE} | awk '{sum+=$3;cnt++}END{print sum/cnt" "sum}';
done
