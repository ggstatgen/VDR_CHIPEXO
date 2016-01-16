#!/bin/bash

#25/7/2014
#This processes a bam files by removing all reads which fall in a bed file you provide.
#I created it to remove the encode black listed regions from my alignments prior to computing autocorrelation and other measures from my files.

#command line 
# bedtools intersect -abam VDR_GM06986_reads.bam -b /net/isi-scratch/giuseppe/indexes/ENCODE_BLACKLIST/hg19-wgEncodeDacMapabilityConsensusExcludable.bed -v > VDR_GM06986_reads_RBL.bam

if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH_BAMS> <BED_FILE> <ID>"
        echo "<PATH_BAM> directory for the bam files to clean up"
	echo "<BED_FILE> full path to the bedfile name to use for cleaning up the bams"
	echo "<ID> indicate reference (eg hg19)"
	echo "NOTE BAM AND BED MUST REFER TO THE SAME REFERENCE (eg HG19)"
        exit
fi

PDATA=$1;
PBED=$2;
PCODE="/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin/bedtools";
PID=$3;

for FILE in ${PDATA}/*.bam;
	do ID=`basename ${FILE} ".bam"`;

        SCRIPT=sub_bed_f_bam_${PID}_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	
	#echo "source activate" >> ${PDATA}/${SCRIPT};
	echo "${PCODE} intersect -abam ${FILE} -b ${PBED} -v > ${PDATA}/${ID}_RBL.bam" >> ${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/sub_bed_f_bam_${PID}_${ID}.err -o ${PDATA}/sub_bed_f_bam_${PID}_${ID}.out -q fgu217.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
