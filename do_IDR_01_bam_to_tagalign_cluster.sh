#!/bin/bash

#from 
#https://sites.google.com/site/anshulkundaje/projects/idr#TOC-IDR-ANALYSIS-ON-SELF-PSEUDOREPLICATES
#INPUT: bam file from bwa, stampy or bowtie
#OUTPUT: tagAlign file for SPP or phantompeakcall or Anshul Kundaje's other stuff

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "<PATH> data path for the bam files"
        exit
fi

PDATA=$1;
PCODE="/home/giuseppe/local/bin";
PBEDTOOLS='/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin'

for FILE in ${PDATA}/*.bam;
        do ID=`basename ${FILE} ".bam"`;

        SCRIPT=bam2tgln_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate' >>${PDATA}/${SCRIPT};


	echo "${PCODE}/samtools view -b -F 1548 -q 30 ${FILE} | ${PBEDTOOLS}/bamToBed -i stdin | awk 'BEGIN{FS=\"\t\";OFS=\"\t\"}{\$4=\"N\"; print \$0}' | gzip -c > ${PDATA}/${ID}.tagAlign.gz"	>>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/bam2tgln_${ID}.err -o ${PDATA}/bam2tgln_${ID}.out  -q medium_jobs.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT};
done

