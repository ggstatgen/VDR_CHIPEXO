#!/bin/bash

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
	echo "NOTE: using Encode Blacklist regions under /net/isi-scratch/giuseppe/indexes/BLACKLIST_ALLTRACKS/b37_anshuldac_duke.bed";
        exit
fi

PDATA=$1;
ALLSITES=$2;
QVAL=$3
PCODE="/net/isi-backup/giuseppe/scripts/P_Various/do_alleleseq_pp_01_filter_events_in_EBL.pl";

for FILE in ${PDATA}/interestingHets*.txt;
        do
	ID=`echo ${FILE} | grep -Po "NA\d{5}"`;
	BASENAME=`basename ${FILE} ".txt"`;
	
	SCRIPT=asb_filter_EBL_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate' >>${PDATA}/${SCRIPT};

	echo "${PCODE} -i=${FILE}" >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/asb_filter_EBL_${ID}.err -o ${PDATA}/asb_filter_EBL_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
