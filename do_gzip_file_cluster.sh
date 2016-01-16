#!/bin/bash
#script to gzip text files in parallel using the cluster

if [ ! $# == 2 ]; then
	echo "Usage: `basename $0` <PATH> <EXTENSION>"
        echo "You must specify 1) data path (e.g. /home/me/files) and 2) file extension (e.g. fastq)"
	exit
fi

PDATA=$1;
EXT=$2;

for FILE in ${PDATA}/*.${EXT};
        #do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
	do ID=`basename ${FILE} ".${EXT}"`;	

	SCRIPT=job_${ID}.sh;
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "gzip ${FILE}" >>${PDATA}/${SCRIPT};
	nice -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
