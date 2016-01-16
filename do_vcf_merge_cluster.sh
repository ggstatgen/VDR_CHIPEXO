#!/bin/bash

#This script uses a high-speed C library, HTSlib, to merge two VCF
#This is different from concatenating them:
#Here I'm interested in putting together the data I got from hapmap. The pipeline was:
#get hapmap rawfiles, CEU and YRI
#hapmap -> vcf (gatk)
#bzip compress and tabix index
#merge YRI+CEU


if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <DATA1> <DATA2> <OUT>"
        echo "DATA1 - directory containing vcf files to merge with files in DATA2"
        echo "DATA2 - directory containing vcf files to merge with files in DATA1"
        echo "OUT - output data path"
        exit
fi

PDATA1=$1;
PDATA2=$2;
POUT=$3;
BASENAME="genotypes_CEUYRI_r27_nr.b36_fwd"; #basename of the output file, change here

PCODE_C="/net/isi-scratch/giuseppe/tools/htslib/vcf merge"; #C version
PCODE_PERL="/net/isi-scratch/giuseppe/tools/vcftools_0.1.10/bin/vcf-merge" #Perl version - slower but maybe better output

for FILE1 in ${PDATA1}/*.vcf.gz;
        do 
        #BASENAME=`basename ${FILE1} ".vcf"`;
	ID1=`echo ${FILE1} | egrep -o "chr[0-9+|A-Z]*"`;
	for FILE2 in ${PDATA2}/*.vcf.gz;
		do
		ID2=`echo ${FILE2} | egrep -o "chr[0-9+|A-Z]*"`;
		if [ ${ID2} = ${ID1} ]; then
			#run merge
			SCRIPT=vcf_merge_${ID1}_${ID2}.sh;
                        echo '#!/bin/bash' >>${POUT}/${SCRIPT};
                        echo '' >>${POUT}/${SCRIPT};
			echo "${PCODE_PERL} ${FILE1} ${FILE2} | bgzip -c > ${POUT}/${BASENAME}_${ID1}.vcf.gz" >>${POUT}/${SCRIPT};
                        nice -5 qsub -e ${POUT}/vcfmerge_${ID1}.err -o ${POUT}/vcfmerge_${ID1}.out -q newnodes.q ${POUT}/${SCRIPT};
			rm ${OUT}/${SCRIPT};
		fi
	done
done
