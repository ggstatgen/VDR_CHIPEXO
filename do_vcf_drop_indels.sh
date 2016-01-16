#!/bin/bash

#script to remove the indel info from a bunch of vcf
#I want to clean up the 1kg vcfs before getting D' LD blocks with PLIN

#method from here
#http://vcftools.sourceforge.net/man_latest.html
#Output a new vcf file from the input vcf file that removes any indel sites
#vcftools --vcf input_file.vcf --remove-indels --recode --recode-INFO-all --out SNPs_only

#To write out the variants that pass through filters use the --recode option. In addition, use --recode-INFO-all to include all data from the INFO fields in the output. By default INFO fields are not written because many filters will alter the variants in a file, rendering the INFO values incorrect.


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "PATH - path for the vcf files"
        exit
fi

PDATA=$1;
#location of vcf-concat executable, change if needed
PCODE="/home/giuseppe/local/bin/vcftools";

for FILE in ${PDATA}/*.vcf;
        do ID=`basename ${FILE} ".vcf"`;

        SCRIPT=vcftools_drop_indels_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate' >>${PDATA}/${SCRIPT};
        echo "${PCODE} --vcf ${FILE} --remove-indels --recode --recode-INFO-all --out ${PDATA}/${ID}_noindels" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/vcftools_drop_indels_${ID}.err -e ${PDATA}/bcftools_drop_indels_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done

