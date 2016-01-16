#!/bin/bash

PCODE="/net/isi-scratch/giuseppe/tools";

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the bam files"
        exit
fi

PDATA=$1;

for FILE in ${PDATA}/*.bam;
	do 
        BASEFILE=`basename ${FILE} ".bam"`; 
	SCRIPT=${BASEFILE}.sh;
 
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	echo "${PCODE}/samtools/samtools view -f4 ${FILE}" >>${PDATA}/${SCRIPT};
       	nice -5 qsub -e ${PDATA}/${BASEFILE}.err -o ${PDATA}/${BASEFILE}_unmapped.sam -q short_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT}; 
done
