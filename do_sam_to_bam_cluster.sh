#!/bin/bash

#use this if you already have bam files

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input sams"
        exit
fi

PDATA=$1;
PCODE='/home/giuseppe/local/bin';

for FILE in ${PDATA}/*.sam;
        do 
        ID=`basename ${FILE} ".sam"`;
	#ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        SCRIPT=samtools_${ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "${PCODE}/samtools view -bS ${FILE} > ${PDATA}/${ID}.bam;" >> ${PDATA}/${SCRIPT};
        echo "${PCODE}/samtools sort ${PDATA}/${ID}.bam ${PDATA}/${ID}.sorted;" >>${PDATA}/${SCRIPT};
	echo "${PCODE}/samtools index ${PDATA}/${ID}.sorted.bam;" >>${PDATA}/${SCRIPT};
	 
	nice -5 qsub -e ${PDATA}/samtools_${ID}.err -o ${PDATA}/samtools_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done

