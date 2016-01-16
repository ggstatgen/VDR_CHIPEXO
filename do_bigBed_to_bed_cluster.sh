#!/bin/bash
#script to conver bigbed to bed using the UCSC script
#bigbed:
#http://genome.ucsc.edu/goldenPath/help/bigBed.html

#/net/isi-scratch/giuseppe/tools/UCSC_tools/bigBedToBed 
#bigBedToBed - Convert from bigBed to ascii bed format.
#usage:
#   bigBedToBed input.bb output.bed
#options:
#   -chrom=chr1 - if set restrict output to given chromosome
#   -start=N - if set, restrict output to only that over start
#   -end=N - if set, restict output to only that under end
#   -maxItems=N - if set, restrict output to first N items
#   -udcDir=/dir/to/cache - place to put cache for remote bigBed/bigWigs


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the *compressed* wig files"
        exit
fi

PDATA=$1;
#location of wigtobigwig executable, change if needed
PCODE="/net/isi-scratch/giuseppe/tools/UCSC_tools";

for FILE in ${PDATA}/*.bb;
	do 
        BASEFILE=`basename ${FILE} ".bb"`; 
        SCRIPT=bb_to_b_${BASEFILE}.sh;
       
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
        echo "${PCODE}/bigBedToBed ${FILE} ${PCHROM} ${PDATA}/${BASEFILE}.bed;" >>${PDATA}/${SCRIPT};
        nice -5 qsub -e ${PDATA}/bb_to_b_${BASEFILE}.err -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT}; 
done
