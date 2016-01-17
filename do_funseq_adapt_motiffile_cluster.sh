#!/bin/bash

#wrapper to run do_funseq_adapt_motiffile.pl on the cluster

#Input: a directory full of files created by do_pscanchip_out_intersect_vdrbvs.pl, in PscanCHIP format, each for a different TF whose PWM intervals around your sequences you've chosen to download

#This script will convert those files into a format suitable for Funseq. Collate the outputs of this and rename "ENCODE.tf.bound.union.bed", then save in the /data dir of Funseq2 (backing up the previous version)

#synopsis
#	print "USAGE: do_funseq_adapt_motiffile.pl -i=<INFILE> -pwm=<ENCODE_PWM_FILE> -id=<ID> -m=<MOTIF_NAME>\n";
#	print "<INFILE> .ris file from PscanChIP or .out file from do_pscanchip_out_intersect_vdrbvs.pl\n";
#	print "<ENCODE_PWM_FILE> text file with Jaspar PWMs in ENCODE format obtained with RSAT\n";
#	print "<ID> string to use for the ID (e.g. Jaspar ID) of the PWM in the output file\n";
#	print "<MOTIF_NAME> string to use for the Motif PWM name in the output file\n";


if [ ! $# == 4 ]; then
	echo "Usage: `basename $0` <VDR_BV_PATH> <RIS_PATH> <ENCODE_PWM_FILE> <EXT> <MINSCORE>"
	echo "<VDR_BV_PATH> full path to the Output_noDBRECUR.vcf file of variants"
	echo "<RIS_PATH> - absolute path for the pscanchip.ris files, one per TF"
	echo "<ENCODE_PWM_FILE> text file with Jaspar PWMs in ENCODE format obtained with RSAT"
	echo "<MINSCORE> pscanchip score threshold for minimun motif score. NA if argument not needed"
	echo "NOTE: this will only work with filenames with the structure Pscanchip_allsites_BVs_<TF>_<ID>_sites.<EXT>";
	echo "so modify the filename accordingly";
	exit
fi

PVCF=$1
PDATA=$2;
PPWM=$3;
MINSCORE=$4;
PCODE="/net/isi-backup/giuseppe/scripts/P_Various/do_funseq_adapt_motiffile.pl";

for FILE in ${PDATA}/Pscanchip*.ris;
	do
	SEED=`echo ${FILE} | grep -Po "BVs_(.*)_sites"`;
	TF=`echo ${SEED} | cut -d _ -f 2`;
	ID=`echo ${SEED} | cut -d _ -f 3`;

	if [ ${MINSCORE} == 'NA' ]; then
		SCRIPT=script_adapt_motiffile_${TF}_${ID}_mscoreNA.sh;
	else
		SCRIPT=script_adapt_motiffile_${TF}_${ID}_mscore${MINSCORE}.sh;
	fi

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	if [ ${MINSCORE} == 'NA' ]; then 
		echo "perl ${PCODE} -i_vcf=${PVCF} -i_ris=${FILE} -pwm=${PPWM} -id=${ID} -m=${TF}" >>${PDATA}/${SCRIPT};
		nice -5 qsub -e ${PDATA}/script_adapt_motiffile_${TF}_${ID}_mscoreNA.err -o ${PDATA}/Pscanchip_funseqadapt_${TF}_${ID}_mscoreNA.out -q medium_jobs.q ${PDATA}/${SCRIPT};
	else
		echo "perl ${PCODE} -i_vcf=${PVCF} -i_ris=${FILE} -pwm=${PPWM} -id=${ID} -m=${TF} -s=${MINSCORE}" >>${PDATA}/${SCRIPT};
		nice -5 qsub -e ${PDATA}/script_adapt_motiffile_${TF}_${ID}_mscore${MINSCORE}.err -o ${PDATA}/Pscanchip_funseqadapt_${TF}_${ID}_mscore${MINSCORE}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
	fi
	rm ${PDATA}/${SCRIPT};  
done
