#!/bin/bash
#script to build tabix index for bzipped files

if [ ! $# == 1 ]; then
	echo "Usage: `basename $0` <PATH>"
        echo "You must specify the data path (e.g. /home/me/files)"
	exit
fi

PDATA=$1;
PCODE="/net/isi-scratch/giuseppe/tools/tabix-0.2.6";
EXT=$2;

for FILE in ${PDATA}/*.vcf.gz;
        #do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;  FILE_ID=`basename ${FILE} ".txt"`;
	do ID=`basename ${FILE} ".vcg.gz"`;	

	SCRIPT=${ID}.sh;
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "${PCODE}/tabix ${FILE}" >>${PDATA}/${SCRIPT};
	nice -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
