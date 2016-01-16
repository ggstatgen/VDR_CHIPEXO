#!/bin/bash

#Script to concatenate VCF files (eg from different chrs)
#from
#http://vcftools.sourceforge.net/perl_module.html#vcf-concat

#syntax
#vcf-concat A.vcf.gz B.vcf.gz C.vcf.gz | gzip -c > out.vcf.gz


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "PATH - path for the vcf.gz files"
        exit
fi

PDATA=$1;
#location of vcf-concat executable, change if needed
PCODE="/home/giuseppe/local/bin/vcf-concat";


SPACE=" ";
for FILE in ${PDATA}/*.vcf.gz;
        do
        INPUT+=${FILE}${SPACE};
done

SCRIPT=vcf-concat.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "${PCODE} ${INPUT} | bgzip -c > ${PDATA}/out.vcf.gz" >> ${PDATA}/${SCRIPT};

nice -5 qsub -e ${PDATA}/vcf-concat.err -o ${PDATA}/vcf-concat.out -q newnodes.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};
