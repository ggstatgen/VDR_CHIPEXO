#!/bin/bash
#create quick bai index for bam alignment, used to visualise bam in IGV

if [ ! $# == 1 ]; then
	echo "Usage: `basename $0` <PATH>"
        echo "PATH - directory containing the bams"
	exit
fi

PDATA=$1;
PCODE="/home/giuseppe/local/bin";
EXT=$2;

for FILE in ${PDATA}/*.bam;
        do 
	#ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;  
        #FILE_ID=`basename ${FILE} ".txt"`;
	ID=`basename ${FILE} ".bam"`;	

	SCRIPT=samtools_idx_${ID}.sh;
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "${PCODE}/samtools index ${FILE}" >>${PDATA}/${SCRIPT};
	nice -5 qsub -e ${PDATA}/smt_idx_${ID}.err -o ${PDATA}/smt_idx_${ID}.out -q fgu217.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
