#!/bin/bash

#command: /net/isi-scratch/giuseppe/tools/bwa-0.6.2/bwa aln -t8 -f [output.sai] [bwa_genome_index] [reads.fq]
#seguito da
#/net/isi-scratch/giuseppe/tools/bwa-0.6.2/bwa samse [bwa_genome_index] [output.sai] [reads.fq] -f [output.sam]

PDATA='/net/isi-scratch/giuseppe/VDR/MY_BAM_STAMPY/';
PCODE='/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin/';

for FILE in ${PDATA}*.bam;
   do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
   ${PCODE}macs14 -t ${FILE} -n ${PDATA}${ID} -g hs --keep-dup=all -w -S;
done
