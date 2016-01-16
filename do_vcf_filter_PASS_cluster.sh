#!/bin/bash
#this calls the perl script clean_vcf.pl
#and shares the workload on the cluster

#USAGE: clean_vcf.pl -i=<INFILE>
#<INFILE> .vcf file to clean


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input .vcf.gz files"
        exit
fi

PDATA=$1;
PCODE="/net/isi-backup/giuseppe/scripts";

for FILE in ${PDATA}/*.vcf.gz;
	do ID=`basename ${FILE} ".vcf.gz"`;
	SCRIPT=do_VCFP_${ID}.sh;

	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "gunzip ${FILE}" >>${PDATA}/${SCRIPT};
	echo "perl ${PCODE}/clean_vcf.pl -i=${ID}.vcf" >>${PDATA}/${SCRIPT};
	echo "gzip ${PDATA}/${ID}.vcf" >>${PDATA}/${SCRIPT};
        nice -5 qsub -v "BASH_ENV=~/.bashrc" -e ${PDATA}/VCF_P_${ID}.err -o ${PDATA}/VCF_P_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
