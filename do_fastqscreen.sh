#!/bin/bash

PCODE="/net/isi-scratch/giuseppe/tools";
if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the fastq.gz files"
        exit
fi

PDATA=$1;

for FILE in ${PDATA}/*.fastq;
        do
        PAIR=`echo ${FILE} | egrep -o "_2_"`;
        if [ "${PAIR}" = '_2_' ];
        then
                continue;
        fi

        BASEFILE=`basename ${FILE} "_1_filtered_sequence.fastq"`;
        
        ${PCODE}/FastQScreen_v0.3.1/fastq_screen --nohits --illumina1_3 --paired --conf=${PDATA}/fastq_screen.conf  ${PDATA}/${BASEFILE}_1_filtered_sequence.fastq ${PDATA}/${BASEFILE}_2_filtered_sequence.fastq;
done
