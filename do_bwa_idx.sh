#!/bin/bash

PBWA="/home/giuseppe/local/bin" #this refers to a symbolic link to the latest bwa. If you don't use this, it will try to use the cgat bwa

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <DATA_DIR>"
        echo "DATA_DIR - full path for the reference genome(s) to index with BWA"
        exit
fi

PDATA=$1;

for FILE in ${PDATA}/*.fa;
        do
        BASEFILE=`basename ${FILE} ".fa"`;
        SCRIPT=${BASEFILE}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        echo "${PBWA}/bwa index -a bwtsw ${FILE}" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/${BASEFILE}.err -o ${PDATA}/${BASEFILE}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
