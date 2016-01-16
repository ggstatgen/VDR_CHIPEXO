#!/bin/bash
#uses samtools to get sam from bams

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input bams"
        exit
fi

PDATA=$1;
PCODE='/net/isi-scratch/giuseppe/tools';

for FILE in ${PDATA}/*.bam;
        do
	#ID=`echo ${FILE} | egrep -o "NA[0-9]*"`;
	ID=`basename ${FILE} ".bam"`;
	
	SCRIPT=bam2sam_${ID}.sh;
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "${PCODE}/samtools/samtools view -h -o ${PDATA}/${ID}.sam ${FILE};" >>${PDATA}/${SCRIPT};
	nice -5 qsub -e ${PDATA}/b2s_${ID}.err -o ${PDATA}/b2s_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
