#!/bin/bash

#run this after running the alleleseq or QTL pipelines and getting LD intervals from SNAP or haploreg.
#This looks if any of the snps in those lists are in disease hits.

#It will run both do_LD_SNP_get_disease_assoc_GWAScatalog.pl and do_LD_SNP_get_disease_assoc_GRASP.pl

#perl do_LD_SNP_get_disease_assoc_GRASP.pl / do_LD_SNP_get_disease_assoc_GWAS.pl
#USAGE: do_LD_SNP_get_disease_assoc.GWAS.pl -i=<LD_SNPS> -source=[snap|haploreg] -p=<SNP_PVAL>
#<LD_SNPS> tab separated output of SNAP/Haploreg for snps proxy search
#[snap|haploreg] source of LD info
#(optional)<SNP_PVAL> significance level for association (default = 5*10-8)

#input is like
#haploreg_0.8_CEU_NA06986.txt

if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <PATH> <EXT> <SOURCE> <PVAL>"
        echo "<PATH> - absolute path for the interestingHets.txt files"
	echo "<EXT> - file extension for the input files"
	echo "<SOURCE> - what are the inputs, one of [snap|haploreg]"
	echo "<PVAL> pval threshold for the resulting hits, if NA default will be used (5*10-8)"
        exit
fi

PDATA=$1;
EXT=$2;
SOURCE=$3;
PVAL=$4;
PCODE="/net/isi-backup/giuseppe/scripts/P_Various";

for FILE in ${PDATA}/*.${EXT};
        do
	#ID=`echo ${FILE} | grep -Po "NA\d{5}"`;
	#ID=$RANDOM;
	BASENAME=`basename ${FILE} ".${EXT}"`;
	ID=`basename ${FILE} ".${EXT}"`;

	if [ $PVAL == 'NA' ]; then 
		SCRIPT=get_da_${ID}.sh;
	else
		SCRIPT=get_da_pthrs${PVAL}_${ID}.sh;
	fi

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	if [ $PVAL == 'NA' ]; then
		echo "${PCODE}/do_LD_SNP_get_disease_assoc_GRASP.pl -i=${FILE} -source=${SOURCE} > ${PDATA}/GRASP_${BASENAME}.txt" >>${PDATA}/${SCRIPT};
		echo "${PCODE}/do_LD_SNP_get_disease_assoc_GWAScatalog.pl  -i=${FILE} -source=${SOURCE} > ${PDATA}/GWAScat_${BASENAME}.txt" >>${PDATA}/${SCRIPT};
	else
		echo "${PCODE}/do_LD_SNP_get_disease_assoc_GRASP.pl -i=${FILE} -source=${SOURCE} -p=${PVAL} > ${PDATA}/GRASP_${BASENAME}.txt" >>${PDATA}/${SCRIPT};
		echo "${PCODE}/do_LD_SNP_get_disease_assoc_GWAScatalog.pl  -i=${FILE} -source=${SOURCE} -p=${PVAL} > ${PDATA}/GWAScat_${BASENAME}.txt" >>${PDATA}/${SCRIPT};
	fi

	if [ $PVAL == 'NA' ]; then
		nice -5 qsub -e ${PDATA}/get_da_${ID}.err -o ${PDATA}/get_da_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
	else
		nice -5 qsub -e ${PDATA}/get_da_pthrs${PVAL}_${ID}.err -o ${PDATA}/get_da_${ID}_qthrs${PVAL}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
	fi
        rm ${PDATA}/${SCRIPT};  
done
