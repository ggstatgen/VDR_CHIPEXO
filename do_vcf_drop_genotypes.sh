#!/bin/bash

#script to remove the genotype info from a vcf
#I created it to slim down the vcfs downloaded from 1kg to be used to create a background for the FIGURE5 analyses

#uses bcftools

#syntax
#About:   VCF/BCF conversion, view, subset and filter VCF/BCF files.
#Usage:   bcftools view [options] <in.vcf.gz> [region1 [...]]
#
#Output options:
#    -G,   --drop-genotypes              drop individual genotype information (after subsetting if -s option set)


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "PATH - path for the vcf.gz files"
        exit
fi

PDATA=$1;
#location of vcf-concat executable, change if needed
PCODE="/home/giuseppe/local/bin/bcftools view -G";

for FILE in ${PDATA}/*.vcf.gz;
        do ID=`basename ${FILE} ".vcf.gz"`;

        SCRIPT=bcftools_drop_gts_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "${PCODE} ${FILE} | gzip -c > ${PDATA}/${ID}_nogt.vcf.gz" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/bcftools_drop_gts_${ID}.err -e ${PDATA}/bcftools_drop_gts_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done

