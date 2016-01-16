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

if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <PATH> <ALLSITES> <QVAL> <MOTIF_FILE>"
        echo "<PATH> - absolute path for the interestingHets.txt files"
	echo "<ALLSITES> [y|n] whether you want to keep in the output ASB events which did not fall in the peak interval specified during Alleleseq analysis";
	echo "<QVAL> qval threshold (if you want to restrict the returned sites, eg 0.05). NA if no threshold";
	echo "<MOTIF_FILE> bed motif file with scores, from PScanChip. Only instances intersecting motifs will be kept. NA if no motif file";
        exit
fi

PDATA=$1;
ALLSITES=$2;
QVAL=$3;
MOTIFS=$4;
PCODE="/net/isi-backup/giuseppe/scripts/P_Various/do_alleleseq_pp_03_get_rsids_from_interestinghets_consensus.pl";

for FILE in ${PDATA}/interestingHets_consensus*.txt;
        do
	#ID=`echo ${FILE} | grep -Po "o\d+"`;
	ID=`basename ${FILE} ".txt"`;


	if [ $QVAL == 'NA' ]; then 
		if [ $MOTIFS == 'NA' ]; then
			SCRIPT=asb_getrsID_${ALLSITES}_${ID}.sh;
		else
			SCRIPT=asb_getrsID_inmotifs_${ALLSITES}_${ID}.sh;
		fi
	else
		if [ $MOTIFS == 'NA' ]; then
			SCRIPT=asb_getrsID_${ALLSITES}_qthrs${QVAL}_${ID}.sh;
		else
			SCRIPT=asb_getrsID_inmotifs_${ALLSITES}_qthrs${QVAL}_${ID}.sh;
		fi
	fi

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	if [ $ALLSITES == 'y' ]; then
		if [ $QVAL == 'NA' ]; then
			if [ $MOTIFS == 'NA' ]; then
				echo "${PCODE} -i=${FILE} -allsites" >>${PDATA}/${SCRIPT};
			else
				echo "${PCODE} -i=${FILE} -allsites -m=${MOTIFS}" >>${PDATA}/${SCRIPT};
			fi
		else
			if [ $MOTIFS == 'NA' ]; then
				echo "${PCODE} -i=${FILE} -allsites -qval=${QVAL} " >>${PDATA}/${SCRIPT};
			else
				echo "${PCODE} -i=${FILE} -allsites -qval=${QVAL}  -m=${MOTIFS}" >>${PDATA}/${SCRIPT};
			fi
		fi
	elif [ $ALLSITES == 'n' ]; then
		if [ $QVAL == 'NA' ]; then
			if  [ $MOTIFS == 'NA' ]; then
				echo "${PCODE} -i=${FILE}" >>${PDATA}/${SCRIPT};
			else
				echo "${PCODE} -i=${FILE} -m=${MOTIFS}" >>${PDATA}/${SCRIPT};
			fi
		else
			if  [ $MOTIFS == 'NA' ]; then
				echo "${PCODE} -i=${FILE} -qval=${QVAL}" >>${PDATA}/${SCRIPT};
			else
				echo "${PCODE} -i=${FILE} -qval=${QVAL}  -m=${MOTIFS}" >>${PDATA}/${SCRIPT};
			fi
		fi
	else
		echo "ERROR: <ALLSITES> variable not recognised. Aborting..";
		rm ${PDATA}/${SCRIPT};
		exit $?
	fi

	if [ $QVAL == 'NA' ]; then
		if  [ $MOTIFS == 'NA' ]; then
			nice -5 qsub -e ${PDATA}/asb_getrsID_${ALLSITES}_${ID}.err -o ${PDATA}/rsID_${ID}_${ALLSITES}.tsv -q medium_jobs.q ${PDATA}/${SCRIPT};
		else
			nice -5 qsub -e ${PDATA}/asb_getrsID_inmotifs_${ALLSITES}_${ID}.err -o ${PDATA}/rsID_${ID}_inmotifs_${ALLSITES}.tsv -q medium_jobs.q ${PDATA}/${SCRIPT};
		fi
	else
		if  [ $MOTIFS == 'NA' ]; then
			nice -5 qsub -e ${PDATA}/asb_getrsID_${ALLSITES}_qthrs${QVAL}_${ID}.err -o ${PDATA}/rsID_${ID}_${ALLSITES}_qthrs${QVAL}.tsv -q medium_jobs.q ${PDATA}/${SCRIPT};
		else
			nice -5 qsub -e ${PDATA}/asb_getrsID_inmotifs_${ALLSITES}_qthrs${QVAL}_${ID}.err -o ${PDATA}/rsID_${ID}_inmotifs_${ALLSITES}_qthrs${QVAL}.tsv -q medium_jobs.q ${PDATA}/${SCRIPT};
		fi
	fi
        rm ${PDATA}/${SCRIPT};  
done
