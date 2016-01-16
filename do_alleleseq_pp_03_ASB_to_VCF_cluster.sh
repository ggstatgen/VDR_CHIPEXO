#!/bin/bash

#run this after you've ran the alleleseq pipeline on multiple single samples
#this returns a VCF 4.1 of the alleleseq (overlap) output.

#So once you have overlap data across samples, this will give you a .vcf version
#you can then do all sorts of arithmetics on it using bedtools (intersect with motifs etc)
#the alleleseq pipeline output , interestinHets.txt, does not contain the SNP ids from 1000g
#this extracts them so you can query dbs like SNAP or HAPLOREG with the output

#USAGE: do_alleleseq_pp_03_ASB_to_VCF.pl -i=<INFILE_ASSOC> -allsites
#<INFILE_ASSOC> file interestingHets_overlaps.txt from do_alleleseq_pp_02..
#(opt)<allsites> flag; if set, retrieve asb sites even when they did not fall in peak interval file.

#NOTE
#HUGE MEMORY REQUIREMENTS
#ONLY RUN ON 217

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <ALLSITES>"
        echo "<PATH> - absolute path for the interestingHets.txt files"
	echo "<ALLSITES> [y|n] whether you want to keep in the output ASB events which did not fall in the peak interval specified during Alleleseq analysis";
        exit
fi

PDATA=$1;
ALLSITES=$2;
PCODE="/net/isi-backup/giuseppe/scripts/P_Various/do_alleleseq_pp_03_ASB_to_VCF.pl";

for FILE in ${PDATA}/interestingHets_consensus*.txt;
        do
	#ID=`echo ${FILE} | grep -Po "o\d+"`;
	ID=`basename ${FILE} ".txt"`;
	SCRIPT=asb2vcf_${ALLSITES}_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	if [ $ALLSITES == 'y' ]; then
		echo "${PCODE} -i=${FILE} -allsites" >>${PDATA}/${SCRIPT};
	elif [ $ALLSITES == 'n' ]; then
		echo "${PCODE} -i=${FILE}" >>${PDATA}/${SCRIPT};
	else
		echo "ERROR: <ALLSITES> variable not recognised. Aborting..";
		rm ${PDATA}/${SCRIPT};
		exit $?
	fi

	nice -5 qsub -e ${PDATA}/asb2vcf_${ALLSITES}_${ID}.err -o ${PDATA}/asb2vcf_${ID}_${ALLSITES}.out -q fgu217.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
