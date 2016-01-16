#!/bin/bash

DATA="/net/isi-scratch/giuseppe/NOTCH_DATA/BAM/";
CODE="/net/isi-scratch/giuseppe/tools/FastQC/";

for FILE in ${DATA}*.bam; 
   do ${CODE}fastqc --outdir=${DATA}FASTQC_REPORTS --threads 6 ${FILE};  
done
