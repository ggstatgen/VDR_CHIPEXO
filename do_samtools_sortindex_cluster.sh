#!/bin/bash

#use this if you already have bam files

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input unsorted bams"
        exit
fi

PDATA=$1;
PCODE='/home/giuseppe/local/bin';

for FILE in ${PDATA}/*.bam;
        do 
        BASENAME=`basename ${FILE} ".bam"`;
	ID=`basename ${FILE} ".bam"`
	#ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        SCRIPT=samtools_sortindex_${ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "${PCODE}/samtools sort ${FILE} ${PDATA}/${BASENAME}.sorted;" >>${PDATA}/${SCRIPT};
	echo "${PCODE}/samtools index ${PDATA}/${BASENAME}.sorted.bam;" >>${PDATA}/${SCRIPT};
	 
	nice -5 qsub -e ${PDATA}/samtools_${ID}.err -o ${PDATA}/samtools_${ID}.out -q fgu217.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done

