#!/bin/bash

#wrapper to run do_REVIEW_get_VDRBV_distrib_around_risPWM.pl on the cluster
#Input: a directory full of files created by do_pscanchip_out_intersect_vdrbvs.pl, in PscanCHIP format, each for a different TF whose PWM intervals around your sequences you've chosen to download
#This script will collect the distances of each VDR-BV or VDR-rBV from the closest motif instance, and plot both a line profile and a histogram

#synopsis
#USAGE: do_REVIEW_get_VDRBV_distrib_around_risPWM.pl -m=<INFILE_PSCANCHIP> -v=<BV|rBV> -t=<THRS> -p=<pdf|svg|jpg|png> (opt)-s=<MINSCORE>
#<INFILE_PSCANCHIP> ris file from PscanChip
#<BV|rBV> if BV, all VDRBV will be used; if rBV, only VDR-rBV will be used
#<THRS> Distance threshold cutoff for plot
#<PLOT> type of R output plot desired
#optional <MINSCORE> lower threshold on score (eg 0.8) (default:none)

if [ ! $# == 5 ]; then
	echo "Usage: `basename $0` <RIS_DIR> <BV|rBV> <THRS> <pdf|svg|jpg|png> <MINSCORE|NA>"
	echo "<RIS_DIR> - absolute path for the pscanchip files (.ris extension), one per TF"
	echo "<BV|rBV> compute distances from VDR-BVs (~40,000) or VDR-rBVs (~350)?"        
	echo "<THRS> distance threshold for thresholded plot: max distance from PWM to show on x-axis of plots"
	echo "<PLOT> one of pdf,svg,jpg,png";
	echo "<MINSCORE|NA> minimum pscanchip score to consider. NA if not required";
	echo "paths to chrom sizes, VDR-BVs, VDR-rBVs bed are HARDCODED. Change if required.";
    exit
fi

PDATA=$1;
PVDR=$2;
PTHRS=$3;
PPLOT=$4;
PSCORE=$5;
PCODE="/net/isi-backup/giuseppe/scripts/P_Various/do_REVIEW_get_VDRBV_distrib_around_risPWM.pl";

for FILE in ${PDATA}/Pscanchip*.ris;
    do
	SEED=`echo ${FILE} | grep -Po "BVs_(.*)_sites"`;
	TF=`echo ${SEED} | cut -d _ -f 2`;
	ID=`echo ${SEED} | cut -d _ -f 3`;
	SCRIPT=script_PWMdist_${TF}_${ID}_VDR-${PVDR}_thrs${PTHRS}_plot${PPLOT}_score${PSCORE}.sh;

    echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
    echo '' >>${PDATA}/${SCRIPT};
	#echo 'source activate'  >>${PDATA}/${SCRIPT};

	if [ "${PSCORE}" = "NA"  ]
	then
		echo "/net/isi-cgat/ifs/apps/apps/perl-5.16.1/bin/perl ${PCODE} -m=${PDATA} -v=${PVDR} -t=${PTHRS} -p=${PPLOT}" >>${PDATA}/${SCRIPT};
	else
		echo "/net/isi-cgat/ifs/apps/apps/perl-5.16.1/bin/perl ${PCODE} -m=${PDATA} -v=${PVDR} -t=${PTHRS} -p=${PPLOT} -s=${PSCORE}" >>${PDATA}/${SCRIPT};
	fi
	
	nice -5 qsub -e ${PDATA}/PWMdist_${TF}_${ID}_VDR-${PVDR}_thrs${PTHRS}_plot${PPLOT}_score${PSCORE}.err -o ${PDATA}/PWMdist_${TF}_${ID}_VDR-${PVDR}_thrs${PTHRS}_plot${PPLOT}_score${PSCORE}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
    rm ${PDATA}/${SCRIPT};  
done
