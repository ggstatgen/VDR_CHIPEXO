#!/bin/bash

#This uses  the HTSeq library to count reads from a bam file for a provided interval file (eg peaks)
#It just counts raw reads, you could then use the output counts to run DEseq or to feed into some normalisation
#This script runs the python script do_bed_htseq_countreads.py in parallel for a lot of bam files
#The do_bed_htseq_countreads.py was written based on the info in http://www-huber.embl.de/users/anders/HTSeq/doc/counting.html
#and http://www-huber.embl.de/users/anders/HTSeq/doc/count.html

#Input1: directory of bam files
#Input2: bed file with common peaks for all the bams

#Ideally you would then put together the outputs in a matrix
#rows: samples
#columns: peaks

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH_BAM> <PATH_BED_FILE>"
        echo "<PATH_BAM> location for the bam alignments"
        echo "<PATH_BED_FILE> location for the bed peak file"
        exit
fi

PDATA=$1;
INTERVAL=$2;
PCODE='/net/isi-backup/giuseppe/scripts';

for FILE in ${PDATA}/*.bam;
        do ID=`basename ${FILE} ".bam"`;

        SCRIPT=HTSeq_count_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate' >> ${PDATA}/${SCRIPT};
	echo "python ${PCODE}/do_bed_htseq_countreads.py --bam ${FILE} --bed ${INTERVAL}" >> ${PDATA}/${SCRIPT};
	#sort -k1,1V
        nice -5 qsub -e ${PDATA}/HTSeq_count_${ID}.err -o ${PDATA}/${ID}_htseq_counts.txt -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done

