#!/bin/bash
#script to run fastqc with custom options on the cluster

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <EXTENSION>"
        echo "You must specify 1) data path (e.g. /home/me/files) and 2) file extension [fastq.gz|fastq|bam]"
        exit
fi

PDATA=$1;
EXT=$2;
#location of fastqc executable, change if needed
PCODE="/net/isi-scratch/giuseppe/tools";
PPERL="/net/isi-cgat/ifs/apps/apps/perl-5.16.1/bin/perl";

mkdir ${PDATA}/FASTQC_REPORTS;

for FILE in ${PDATA}/*.${EXT};
	do 
        ID=`basename ${FILE} ".${EXT}"`;
	#ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
        SCRIPT=fastqc_${ID}.sh;
             
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
        #removed  --kmers 10, insert if needed
        echo "${PPERL} ${PCODE}/FastQC/fastqc --outdir=${PDATA}/FASTQC_REPORTS ${FILE}" >>${PDATA}/${SCRIPT};
        nice -5 qsub -e ${PDATA}/fastqc_${ID}.err -o ${PDATA}/fastqc_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
