#!/bin/bash

#stampy
#[pathtostampy]stampy.py -g [pathtoreference] -h [pathtoreference] -t8 --bamkeepgoodreads -M [pathtodata]

#ulteriori opzioni
# --readgroup=ID:id,tag:value,...    Set read-group tags (ID,SM,LB,DS,PU,PI,CN,DT,PL)  (SAM format)
# --solexa, --solexaold, --sanger    Solexa read qualities (@-based); pre-v1.3 Solexa; and Sanger (!-based, default)
# --bwamark                          Include/mark BWA-mapped reads with XP:Z:BWA tag (produces more output lines)
# --sensitive

#example read groups
#@RG	ID:WTCHG_46599_02	LB:812/12_MPX	SM:ANDREW	PL:ILLUMINA
#@RG	ID:WTCHG_46599_04	LB:812/12_MPX	SM:JEAN	PL:ILLUMINA


#PDATA='/net/isi-scratch/giulio/Data_from_138.37.175.2/Ped_results';
PDATA='/net/isi-scratch/giuseppe/VDR/RAW/POOLED/02_BAM_BWA';
#PDATA='/net/isi-scratch/giuseppe/VDR/MY_FASTQ_CUTADAPT/BAM_BWA';
PCODE='/net/isi-scratch/giuseppe/tools';
#PIDX='/net/isi-scratch/giuseppe/indexes/bwa/mm9';
PIDX='/net/isi-scratch/giuseppe/indexes/stampy/hg19_stampy';
POUT='/net/isi-scratch/giuseppe/VDR/RAW/POOLED/02_BAM_BWA/STAMPY_OUT';
#POUT='/net/isi-scratch/giuseppe/VDR/MY_FASTQ_CUTADAPT_6_15/BAM_STAMPY';
PBWA="/home/giuseppe/local/bin"


for FILE in ${PDATA}/*.bam;
   do ID=`basename ${FILE} ".bam"`;

   #nice -5 nohup ${PCODE}/stampy-1.0.21/stampy.py --solexa --bwa=/net/isi-scratch/giuseppe/tools/bwa-0.6.2/bwa  -g ${PIDX} -h ${PIDX} -t 15 --bamkeepgoodreads -M ${FILE} > ${POUT}/${ID}.sam; 
   nice -5 nohup ${PCODE}/stampy-1.0.21/stampy.py --bwa=${PBWA}/bwa -g ${PIDX} -h ${PIDX} -t 20 --bamkeepgoodreads -M ${FILE} > ${POUT}/${ID}.sam;

   #now convert sam to bam, order bam, and create index
   #${PCODE}/samtools/samtools view -bS ${POUT}/${ID}.sam > ${POUT}/${ID}.bam;
   #${PCODE}/samtools/samtools sort ${POUT}/${ID}.bam ${POUT}/${ID}.sorted;
   #${PCODE}/samtools/samtools index ${POUT}/${ID}.sorted.bam;
done
