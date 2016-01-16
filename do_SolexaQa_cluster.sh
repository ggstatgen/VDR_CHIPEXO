#!/bin/bash
#uses the genetrack peakcaller to find peaks in .bed files obtained by converting bam files using bedtools

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input fastq files"
        exit
fi

PDATA=$1;
PCODE="/net/isi-scratch/giuseppe/scripts";

for FILE in ${PDATA}/*.fastq;
        do 
        #ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        ID=`basename ${FILE} ".bed"`;

        SCRIPT=${ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "${PCODE}/SolexaQA.pl ${FILE} -v -m -d=${PDATA}/solexaQA" >>${PDATA}/${SCRIPT};
        nice -5 qsub -e ${PDATA}/solexaQA/${ID}.err -o ${PDATA}/solexaQA/${ID}.out -v "BASH_ENV=~/.bashrc" -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
