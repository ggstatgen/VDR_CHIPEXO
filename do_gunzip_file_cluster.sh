#!/bin/bash
#script to gunzip text files in parallel using the cluster

if [ ! $# == 1 ]; then
	echo "Usage: `basename $0` <PATH>"
        echo "You must specify the data path (e.g. /home/me/files)"
	exit
fi

PDATA=$1;

for FILE in ${PDATA}/*.gz;
        #do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
	do ID=`basename ${FILE} ".gz"`;
	
	SCRIPT=gunzip_${ID}.sh;
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "gunzip ${FILE}" >>${PDATA}/${SCRIPT};
	nice -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
