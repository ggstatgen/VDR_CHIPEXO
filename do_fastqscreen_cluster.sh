#!/bin/bash

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "<PATH> data path for the fastq files"
        exit
fi

PDATA=$1;
#PCODE='/net/isi-scratch/giuseppe/tools/FastQScreen_v0.3.1/';
PCODE="/net/isi-scratch/giuseppe/tools/fastq_screen_v0.4.2/";
PPERL="/net/isi-cgat/ifs/apps/apps/perl-5.16.1/bin/perl";

for FILE in ${PDATA}/*.fastq;
	do
	ID=`basename ${FILE} ".fastq"`;
		
	SCRIPT=fqscreen_${ID}.sh;
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "${PPERL} ${PCODE}/fastq_screen --conf=${PDATA}/fastq_screen.conf --threads 5 --nohits ${FILE};" >>${PDATA}/${SCRIPT};
	nice -5 qsub -e ${PDATA}/fqscreen_${ID}.err -o ${PDATA}/fqscreen_${ID}.out -v "BASH_ENV=~/.bashrc" -q fgu217.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
