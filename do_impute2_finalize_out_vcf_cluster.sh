#!/bin/bash

#Runs do_impute2_finalize_out_vcf.pl on the cluster
#qctool produces vcf with 
#1 no sample nanmes in the header
#2 no chromosome name

#this script gets the sample names from a textfile (one line per sample name, no spaces) and the chr name from each vcf file name

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "<PATH> directory containing the vcf produced by IMPUTE2-qctool"
        exit
fi

PDATA=$1;
#location of macs2 executable, change if needed
PCODE="/net/isi-backup/giuseppe/scripts";
PSAMPLES="/net/isi-scratch/giuseppe/VDR/VARIANTS/IMPUTATION/CEUYRI_30_samplenames.txt"


for FILE in ${PDATA}/*.vcf;
	do ID=`basename ${FILE} ".vcf"`;

        SCRIPT=vcf_finalize_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	echo "perl ${PCODE}/do_impute2_finalize_out_vcf.pl -i ${FILE} -sample ${PSAMPLES}" >>${PDATA}/${SCRIPT};
        
	nice -5 qsub -e ${PDATA}/vcf_finalize_${ID}.err -o ${PDATA}/vcf_finalize_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
