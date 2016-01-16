#!/bin/bash

#This extracts a bam from a sra file using their toolkit
#http://eutils.ncbi.nih.gov/Traces/sra/?view=software

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "You must specify the data path for the sra files (e.g. /home/me/files)"
        exit
fi

PDATA=$1;
PCODE="/net/isi-scratch/giuseppe/tools";


for FILE in ${PDATA}/*.sra;
	do ID=`basename ${FILE} ".sra"`;

        SCRIPT=SRA_samdump_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	
	echo "${PCODE}/sratoolkit.2.3.4-2-centos_linux64/bin/sam-dump --unaligned ${FILE} | ${PCODE}/samtools/samtools view -Sb -o ${PDATA}/${ID}.bam -" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/SRA_samdump_${ID}.err -o ${PDATA}/SRA_samdump_${ID}.out -q fgu217.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
