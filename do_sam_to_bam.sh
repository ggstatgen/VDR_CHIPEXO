#!/bin/bash

#command: /net/isi-scratch/giuseppe/tools/bwa-0.6.2/bwa aln -t8 -f [output.sai] [bwa_genome_index] [reads.fq]
#seguito da
#/net/isi-scratch/giuseppe/tools/bwa-0.6.2/bwa samse [bwa_genome_index] [output.sai] [reads.fq] -f [output.sam]

PDATA='/net/isi-scratch/giuseppe/VDR/MY_FASTQ/01/';
PCODE='/net/isi-scratch/giuseppe/tools/';
PIDX='/net/isi-scratch/giuseppe/Hsap_ref/ucsc_hg19/';


for FILE in ${PDATA}*.sam;
   do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

   ${PCODE}samtools/samtools view -bS ${FILE} > ${PDATA}${ID}.bam;
   ${PCODE}samtools/samtools sort ${PDATA}${ID}.bam ${PDATA}${ID}.sorted;
   ${PCODE}samtools/samtools index ${PDATA}${ID}.sorted.bam;
done
