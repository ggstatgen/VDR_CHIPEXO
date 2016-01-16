#!/bin/bash

#7/aug/2013
#remove X and Y mapped reads from bam

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input bam files"
        exit
fi

PDATA=$1;
for FILE in ${PDATA}/*.bam;
	do
	nohup nice -19 samtools view -h ${FILE} | grep -v \"chrX\" | grep -v \"chrY\" | samtools view -bS - >& ${PDATA}/${BASEFILE}_nosex.bam 2> ${PDATA}/${BASEFILE}.err &
done




#        do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
#
#	BASEFILE=`basename ${FILE} ".bam"`;
#	
#	SCRIPT=bam_filter_${ID}.sh;
#	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
#        echo '' >>${PDATA}/${SCRIPT};
#	echo "samtools view -h ${FILE} | grep -v \"chrX\" | grep -v \"chrY\" | samtools view -bS - > ${PDATA}/${BASEFILE}_nosex.bam" >>${PDATA}/${SCRIPT};
#	nice -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
#	rm ${PDATA}/${SCRIPT}; 
#done
