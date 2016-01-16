#!/bin/bash

#13/05/2014
#obtains fasta files from fastq. Fastq must be non compressed

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input fastq files"
        exit
fi

PDATA=$1;

for FILE in ${PDATA}/*.fastq;
        do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

	BASEFILE=`basename ${FILE} ".fastq"`;
	
	SCRIPT=fq_to_fa_${ID}.sh;
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "cat ${FILE} | perl -e '\$i=0; while(<>){if(/^\@/ && \$i==0){s/^\@/\>/; print;} elsif(\$i==1){print; \$i=-3 } \$i++;}' > ${PDATA}/${BASEFILE}.fasta" >>${PDATA}/${SCRIPT};
	nice -5 qsub -e ${PDATA}/fq_to_fa_${ID}.err -o ${PDATA}/fq_to_fa_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
