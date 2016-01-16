#!/bin/bash

#use this if you already have bam files

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input bams"
        exit
fi

PDATA=$1;
PCODE='/home/giuseppe/local/bin';

for FILE in ${PDATA}/*.bam;
        do 
        ID=`basename ${FILE} ".bam"`;

        SCRIPT=st_idxst_${ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "${PCODE}/samtools sort -m 5000000000 ${FILE} ${PDATA}/${ID}.sorted;" >>${PDATA}/${SCRIPT};
        echo "${PCODE}/samtools index ${PDATA}/${ID}.sorted.bam;" >>${PDATA}/${SCRIPT};
        
	nice -5 qsub -e ${PDATA}/st_idxst_${ID}.err -o ${PDATA}/st_idxst_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done

