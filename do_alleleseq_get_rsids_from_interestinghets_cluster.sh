#!/bin/bash

#run this after you've ran the alleleseq pipeline on multiple single samples
#the alleleseq pipeline output , interestinHets.txt, does not contain the SNP ids from 1000g
#this extracts them so you can query dbs like SNAP or HAPLOREG with the output

#USAGE: do_alleleseq_get_rsids_from_interestinghets.pl -i=<INFILE_ASSOC> -id=<SAMPLE_ID> -allsites -qval=<QVAL> -m=<INFILE_MOTIFS>
#<INFILE_ASSOC> file interestingHets.txt from AlleleSeq
#<SAMPLE_ID> sample id in format NAxxxxx for the data
#(opt)<allsites> flag; if set, retrieve asb sites even when they did not fall in peak interval file.
#(opt)<QVAL> qval threshold (if you want to restrict the returned sites, eg 0.05)
#(opt)<INFILE_MOTIFS> file bed containing motif instances. Only HET associations intersecting these will be kept.

#input is like
#interestingHets_NA06986.txt


#MODIFY IF YOU WANT TO INCLUDE MOTIF OPTION

if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <ALLSITES> <QVAL>"
        echo "<PATH> - absolute path for the interestingHets.txt files"
	echo "<ALLSITES> [y|n] whether you want to keep in the output ASB events which did not fall in the peak interval specified during Alleleseq analysis";
	echo "<QVAL> qval threshold (if you want to restrict the returned sites, eg 0.05). NA if no threshold";
	echo "NOTE: using variants under /net/isi-scratch/giuseppe/VDR/VARIANTS/1000g_hg19_b37/d_SPLIT_SINGLE_b37/";
	echo "CHANGE THIS IN THE ORIGINAL SCRIPT if you use other variants (e.g. imputed ones)";
        exit
fi

PDATA=$1;
ALLSITES=$2;
QVAL=$3
PCODE="/net/isi-backup/giuseppe/scripts/P_Various/do_alleleseq_get_rsids_from_interestinghets.pl";

for FILE in ${PDATA}/interestingHets*.txt;
        do
	ID=`echo ${FILE} | grep -Po "NA\d{5}"`;
	BASENAME=`basename ${FILE} ".txt"`;


	if [ $QVAL == 'NA' ]; then 
		SCRIPT=asb_getrsID_allsites_${ALLSITES}_${ID}.sh;
	else
		SCRIPT=asb_getrsID_allsites_${ALLSITES}_qthrs${QVAL}_${ID}.sh;
	fi

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	if [ $ALLSITES == 'y' ]; then
		if [ $QVAL == 'NA' ]; then
			echo "${PCODE} -i=${FILE} -id=${ID} -allsites" >>${PDATA}/${SCRIPT};
		else
			echo "${PCODE} -i=${FILE} -id=${ID} -allsites -qval=${QVAL} " >>${PDATA}/${SCRIPT};
		fi
	elif [ $ALLSITES == 'n' ]; then
		if [ $QVAL == 'NA' ]; then
			echo "${PCODE} -i=${FILE} -id=${ID}" >>${PDATA}/${SCRIPT};
		else
			echo "${PCODE} -i=${FILE} -id=${ID} -qval=${QVAL}" >>${PDATA}/${SCRIPT};
		fi
	else
		echo "ERROR: <ALLSITES> variable not recognised. Aborting..";
		rm ${PDATA}/${SCRIPT};
		exit $?
	fi

	if [ $QVAL == 'NA' ]; then
		nice -5 qsub -e ${PDATA}/asb_getrsID_allsites_${ALLSITES}_${ID}.err -o ${PDATA}/rsID_${BASENAME}_${ALLSITES}.tsv -q medium_jobs.q ${PDATA}/${SCRIPT};
	else
		nice -5 qsub -e ${PDATA}/asb_getrsID_allsites_${ALLSITES}_qthrs${QVAL}_${ID}.err -o ${PDATA}/rsID_${BASENAME}_${ALLSITES}_qthrs${QVAL}.tsv -q medium_jobs.q ${PDATA}/${SCRIPT};
	fi
        rm ${PDATA}/${SCRIPT};  
done
