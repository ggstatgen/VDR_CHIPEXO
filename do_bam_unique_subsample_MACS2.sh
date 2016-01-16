#!/bin/bash

#I want to get peak positions across all samples and look at peak positions where there is total depletion vs peak
#This follows the suggestion of the MACS2 creator Tao Liu:
#1 remove dups from reads using MACS2 - output a bed
#2 randomly subsample lines from this bed
#3 call peaks on these 
#4 intersect peaks across all samples if the peak appears at least N times

#This uses randomlines.py <input.bed> <samples> <output.bed>

#this follows these instructions
#https://groups.google.com/forum/#!msg/macs-announcement/VkIeAn29ZVw/7zBRuvc2EaAJ
#An alternative way is to first get unique locations and then subsample in separate steps. I have a script in MACS2 package called "filterdup" which can convert any format supported by MACS to BED and filter out the duplicate locations. Use it like:
#$ filterdup -t input.bam -o input_unique.bed -g hs --keep-dup 1
#You need to install MACS2 in order to use it. <https://github.com/taoliu/MACS>
#Then use any script which can randomly pick lines from a plain text file, for example the one in Galaxy 
#https://bitbucket.org/galaxy/galaxy-dist/src/949e4f5fa03a/tools/filters/randomlines.py
#$ randomlines.py input_unique.bed 345678 > subsampled_unique.bed
#You will get 345,678 unique reads.

if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <SAMPLE> <QVAL>"
        echo "<PATH> path for the bam files"
	echo "<SAMPLE> size of the final sample of aligned reads in bp (e.g. 3000000)"
	echo "<QVAL> qvalue for peak calling (e.g. 0.05)";
        exit
fi

PDATA=$1;
SAMPLE=$2;
QVAL=$3;
#location of macs2 executable, change if needed
PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";
PSCRIPT="/net/isi-backup/giuseppe/scripts";

POUT="${PDATA}/d_MACS2_sample_${SAMPLE}_q${QVAL}";
mkdir ${POUT};
P_BED="${POUT}/d_bed";
P_NPK="${POUT}/d_npk";
P_XLS="${POUT}/d_xls";
mkdir ${P_BED};
mkdir ${P_NPK};
mkdir ${P_XLS};

for FILE in ${PDATA}/*.bam;
        do ID=`basename ${FILE} ".bam"`;

        SCRIPT=subsample_bam_MACS2_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};	
	
	echo "source activate" >> ${PDATA}/${SCRIPT};
	echo "macs2 filterdup -i ${FILE} -f BAM -g 2540757438 --keep-dup 1 -o ${POUT}/${ID}_unique.bed" >> ${PDATA}/${SCRIPT};
	#randomly subsample
	echo "${PSCRIPT}/randomlines.py ${POUT}/${ID}_unique.bed ${SAMPLE} ${POUT}/${ID}_unique_sub${SAMPLE}.bed" >> ${PDATA}/${SCRIPT};
	#which ones were smaller than the desired sample size?
	echo "grep "asked to select more lines" ${POUT}/*.err > ${POUT}/FAILED.txt" >> ${PDATA}/${SCRIPT};
	echo "rm ${POUT}/${ID}_unique.bed" >>${PDATA}/${SCRIPT};
	#peak calling on the subsampled
	echo "${PCODE}/macs2 callpeak -f BED -t ${POUT}/${ID}_unique_sub${SAMPLE}.bed  -n ${POUT}/MACS2_${ID}_unique_sub${SAMPLE} -g 2540757438 --keep-dup=1 --nomodel --extsize 26 --call-summits" >>${PDATA}/${SCRIPT};
	#housekeeping
	echo "mv ${POUT}/MACS2_${ID}_unique_sub${SAMPLE}*.xls ${P_XLS}" >>${PDATA}/${SCRIPT};
	echo "mv ${POUT}/MACS2_${ID}_unique_sub${SAMPLE}*.narrowPeak ${P_NPK}" >>${PDATA}/${SCRIPT};
	echo "mv ${POUT}/MACS2_${ID}_unique_sub${SAMPLE}*summit* ${P_BED}" >>${PDATA}/${SCRIPT};
	echo "gzip ${POUT}/${ID}_unique_sub${SAMPLE}.bed" >>${PDATA}/${SCRIPT};	

        nice -5 qsub -e ${POUT}/subsample_bam_${ID}.err -o ${POUT}/subsample_bam_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT};
done     
