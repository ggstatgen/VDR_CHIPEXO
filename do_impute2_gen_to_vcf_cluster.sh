#!/bin/bash

#25/4/2014
#this script is used to run qctools
#http://www.well.ox.ac.uk/~gav/qctool/#tutorial
#to

#1 filter poor variants from .gen IMPUTE2 outputs
#2 obtain a .vcf file out of the .gen file returned by IMPUTE2

#the program subsets the snps based on some QC parameters
#the most important are INFO (normally < 0.4 should be discarded) and HWE p < 1e-6
#NOTE I DON'T KNOW WHAT HWE is for, so I don't use it currently.

#see mailing list OXSTATGEN for thread "Post-imputation QC"

#USAGE
#subset snps
#$ qctool -g example.bgen -og subsetted.gen -snp-missing-rate 0.05 -maf 0.01 1 -info 0.9 1 -hwe 20
#then convert to vcf
#$ qctool -g subsetted.gen -og subsetted.vcf

#FULL TUTORIAL
#http://www.well.ox.ac.uk/~gav/qctool/#tutorial

if [ ! $# == 1 ]; then
        #echo "Usage: `basename $0` <PATH> <MIN_INFO> <HWE>"
	echo "Usage: `basename $0` <PATH>"
        echo "<PATH> data path for the .gen files (IMPUTE2 output)"
	#echo "<MIN_INFO> min threshold for info field (0.4 recommended)"
	#echo "<HWE> don't know what this is, check"
        exit
fi

PDATA=$1;
#INFO=$2;
#HWE=$3
PCODE="/net/isi-scratch/giuseppe/tools/qctool_v1.3-linux-x86_64/";

for FILE in ${PDATA}/*.gen;
	do ID=`basename ${FILE} ".gen"`;

        SCRIPT=qctool_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#1 filter
	#removed -hwe 1E-6
	echo "${PCODE}/qctool -g ${FILE} -og ${PDATA}/${ID}_subset.gen -info 0.4 1" >> ${PDATA}/${SCRIPT};

	#eg qctool -g IMPUTE2_chr13.out -og IMPUTE2_chr13.vcf
	echo "${PCODE}/qctool -g ${PDATA}/${ID}_subset.gen  -og ${PDATA}/${ID}.vcf" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/qctool_${ID}.err -o ${PDATA}/qctool_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
