#!/bin/bash


#first call bwa using g1k index
/net/isi-scratch/giuseppe/tools/bwa-0.6.2/bwa aln -t5 -f /net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301.sai /net/isi-scratch/giuseppe/Hsap_ref/g1k_v37/human_g1k_v37.fasta /net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301.fq
#get sam
/net/isi-scratch/giuseppe/tools/bwa-0.6.2/bwa samse /net/isi-scratch/giuseppe/Hsap_ref/g1k_v37/human_g1k_v37.fasta /net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301.sai /net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301.fq -f /net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301_g1k.sam

#remove temp sai file
rm *.sai

#then call bwa using ucsc_19 index
/net/isi-scratch/giuseppe/tools/bwa-0.6.2/bwa aln -t5 -f /net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301.sai /net/isi-scratch/giuseppe/Hsap_ref/ucsc_hg19/hg19.fa /net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301.fq
#get sam
/net/isi-scratch/giuseppe/tools/bwa-0.6.2/bwa samse /net/isi-scratch/giuseppe/Hsap_ref/ucsc_hg19/hg19.fa  /net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301.sai /net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301.fq -f /net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301_ucsc.sam


#create bam
java -Xmx4g -Djava.io.tmpdir=/net/isi-scratch/giuseppe/dump \
-jar /net/isi-scratch/giuseppe/tools/picard-tools-1.79/SortSam.jar \
SO=coordinate \
INPUT=/net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301_g1k.sam \
OUTPUT=/net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301_g1k.bam \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=true

java -Xmx4g -Djava.io.tmpdir=/net/isi-scratch/giuseppe/dump \
-jar /net/isi-scratch/giuseppe/tools/picard-tools-1.79/SortSam.jar \
SO=coordinate \
INPUT=/net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301_ucsc.sam \
OUTPUT=/net/isi-scratch/giuseppe/VDR/MY_FASTQ/VDR_GM06986_Peconic20301_ucsc.bam \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=true

#then call stampy
#/net/isi-scratch/giuseppe/tools/stampy-1.0.20/stampy.py -g /net/isi-scratch/giuseppe/VDR/SVR11-P001_hg18/i_My_FastQ/hg19_capital -h /net/isi-scratch/giuseppe/VDR/SVR11-P001_hg18/i_My_FastQ/hg19_capital -M /net/isi-scratch/giuseppe/VDR/SVR11-P001_hg18/i_My_FastQ/vitaminDreceptor_sc-1008_LCL_GM06986_-_-_XO111_SVR11-P001_Peconic20301hg18.fq -v 2
