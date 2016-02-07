#!/bin/bash
#uses the bedtools utility bamtobed to convert several bams to bed.
#these are needed as the input of the genetrack peakcaller.

#Tool:    bedtools bamtobed (aka bamToBed)
#Version: v2.17.0
#Summary: Converts BAM alignments to BED6 or BEDPE format.
#
#Usage:   bedtools bamtobed [OPTIONS] -i <bam> 
#
#Options: 
#	-bedpe	Write BEDPE format.
#		- Requires BAM to be grouped or sorted by query.
#
#	-bed12	Write "blocked" BED format (aka "BED12").
#
#		http://genome-test.cse.ucsc.edu/FAQ/FAQformat#format1
#
#	-split	Report "split" BAM alignments as separate BED entries.
#
#	-ed	Use BAM edit distance (NM tag) for BED score.
#		- Default for BED is to use mapping quality.
#		- Default for BEDPE is to use the minimum of
#		  the two mapping qualities for the pair.
#		- When -ed is used with -bedpe, the total edit
#		  distance from the two mates is reported.
#
#	-tag	Use other NUMERIC BAM alignment tag for BED score.
#		- Default for BED is to use mapping quality.
#		  Disallowed with BEDPE output.
#
#	-color	An R,G,B string for the color used with BED12 format.
#		Default is (255,0,0).
#
#	-cigar	Add the CIGAR string to the BED entry as a 7th column.


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input bams"
        exit
fi

PDATA=$1;
PCODE="/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/";

for FILE in ${PDATA}/*.bam;
        do ID=`basename ${FILE} ".bam"`;
	
	SCRIPT=bam2bed_${ID}.sh;
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate'  >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "${PCODE}/bamToBed -i ${FILE} > ${PDATA}/${ID}.bed;" >>${PDATA}/${SCRIPT};
	nice -5 qsub -e ${PDATA}/bam2bed_${ID}.err -o ${PDATA}/bam2bed_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
