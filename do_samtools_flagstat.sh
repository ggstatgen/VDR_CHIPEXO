#!/bin/bash


# gather samtools stats using flagstat and idxstat

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "You must specify the data path for the bam files (e.g. /home/me/files)"
        exit
fi

PDATA=$1;
#location of fastqc executable, change if needed
PCODE='/home/giuseppe/local/bin';


for FILE in ${PDATA}/*.bam;
	do
	ID=`basename ${FILE} ".bam"`;
	echo ""
	echo "$ID: "
        ${PCODE}/samtools flagstat ${FILE};
	#${PCODE}/samtools/samtools idxstats ${FILE} | awk 'BEGIN {a=0;b=0} {a += $3; b+=$4 } END{print a " mapped " b " unmapped "}';
done
