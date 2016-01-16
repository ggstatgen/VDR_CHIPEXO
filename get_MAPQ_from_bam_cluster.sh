#!/bin/bash
#open a bam alignment, get its mapq values, save in textfile to be opened by matlab


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input bam files"
        exit
fi

PDATA=$1;
PCODE="/net/isi-scratch/giuseppe/tools";

for FILE in ${PDATA}/*.bam;
        do 
        #ID=`basename ${FILE} ".bam"`;
	ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        SCRIPT=${ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "${PCODE}/samtools/samtools view ${FILE} | cut -f 5 > ${PDATA}/${ID}_MAPQ.dat;" >>${PDATA}/${SCRIPT};
        nice -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.out -v "BASH_ENV=~/.bashrc" -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
