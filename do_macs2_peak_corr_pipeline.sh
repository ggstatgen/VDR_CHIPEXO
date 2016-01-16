#!/bin/bash

#uses macs2 and UCSC tools to obtain correlations between pairs of bigwig files
#1 call peaks with macs2
#2 call signal bdg tracks with macs2 bdgcmp (FE and logLR only)
#3 converts bdg signal tracks into bw signal tracks
#4 compare each bw with each other and get matrix of correlations ready for R

#input: aligned bam files
#output pdf plot of correlation in R

if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <PATH> <CHROM_INFO> <MACSDUP> <MACS_Q>"
        echo "PATH - path for the bam files (e.g. /home/me/files)"
	echo "CHROM_INFO - file of chromosome sizes"
        echo "MACSDUP - one of [all|auto|1]"
        echo "MACS_Q - eg 0.001 (3) or 0.00001 (5)"
        exit
fi

PDATA=$1;
PCHROM=$2;
MACSDUP=$3; #MACSDUP=all;
MACSQ=$4; #MACSQ=0.00001;

PCODEM="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";
PCODE="/net/isi-scratch/giuseppe/tools";
PSCRIPT="/net/isi-backup/giuseppe/scripts";
MAPPABILITY="2540757438"; #hs, 40bp 

#create output dir
POUT=${PDATA}/d_correlation_${MACSDUP}_q${MACSQ};
mkdir ${POUT};

for FILE in ${PDATA}/*.bam;
	do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        SCRIPT=ppl_corr_1_${ID}.sh;

        echo '#!/bin/bash' >>${POUT}/${SCRIPT};
        echo '' >>${POUT}/${SCRIPT};
	#----------------------------------------------------------------------------------------------------
	#1. get peaks
        #----------------------------------------------------------------------------------------------------
	#change required options here:
        echo "${PCODEM}/macs2 callpeak -t ${FILE} -n ${POUT}/${ID} -g ${MAPPABILITY} --keep-dup=${MACSDUP}  -B --SPMR -q ${MACSQ} --nomodel --extsize 26" >> ${POUT}/${SCRIPT};
        #These are MACS2 outputs as of 17/7/2013
        #$ID_control_lambda.bdg
        #$ID_peaks.bed
        #$ID_peaks.narrowPeak
        #$ID_peaks.xls
        #$ID_summits.bed
        #$ID_treat_pileup.bdg
	#delete all but the bdg
	echo "rm ${POUT}/${ID}_peaks.narrowPeak" >> ${POUT}/${SCRIPT};
	echo "rm ${POUT}/${ID}_peaks.bed" >> ${POUT}/${SCRIPT};
	echo "rm ${POUT}/${ID}_peaks.xls" >> ${POUT}/${SCRIPT};
	echo "rm ${POUT}/${ID}_summits.bed" >> ${POUT}/${SCRIPT};
	#---------------------------------------------------------------------------------------------------
	#2. get bdg signal tracks (FE, logLR, ppois)
	#---------------------------------------------------------------------------------------------------
        #change required options here:
        echo "${PCODEM}/macs2 bdgcmp -t ${POUT}/${ID}_treat_pileup.bdg -c ${POUT}/${ID}_control_lambda.bdg -o ${POUT}/${ID}_ppois.bdg" >>${POUT}/${SCRIPT};
        echo "${PCODEM}/macs2 bdgcmp -t ${POUT}/${ID}_treat_pileup.bdg -c ${POUT}/${ID}_control_lambda.bdg -o ${POUT}/${ID}_FE.bdg -m FE" >>${POUT}/${SCRIPT};
        echo "${PCODEM}/macs2 bdgcmp -t ${POUT}/${ID}_treat_pileup.bdg -c ${POUT}/${ID}_control_lambda.bdg -o ${POUT}/${ID}_logLR.bdg -m logLR -p 0.00001" >>${POUT}/${SCRIPT};

        #remove control data
        echo "rm ${POUT}/${ID}_control_lambda.bdg" >>${POUT}/${SCRIPT};

        nice -5 qsub -e ${POUT}/${ID}_macs.log -q newnodes.q ${POUT}/${SCRIPT};
        rm ${POUT}/${SCRIPT};
done

#you have 4 bdg files for each initial bam
#check that the former set of jobs has finished
#sleep if there are still jobs executing
JOB_STATUS=$(qstat);
if [[ $JOB_STATUS =~ ppl_corr_1 ]] ; 
then 
	COMPLETE=0;
else
	COMPLETE=1; 
fi

while [ "${COMPLETE}" = "0" ]
do
	echo "not all ppl_corr_1 jobs complete. Waiting 60 seconds.."
	sleep 60;

	JOB_STATUS=$(qstat);
	if [[ $JOB_STATUS =~ ppl_corr_1 ]] ; 
	then 
		COMPLETE=0; 
	else
		COMPLETE=1;
	fi
done
echo "All jobs completed, continuing to stage 3..";
echo " ";

#check stderr dumps for errors
#ERR=`egrep '(error|warning|err|warn|fail)' ${PDATA}/*.err`;
#if [[ "$ERR" -ne "\n"  ]] ;
#then
#        echo "Errors in stage 1 .err files. Aborting..";
#	exit -1;
#else
#        rm ${PDATA}/*.err;
#fi

#----------------------------------------------------------------------------------------------------------
#3. bdg to bw
#----------------------------------------------------------------------------------------------------------
for FILE in ${POUT}/*.bdg;
        do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
        BASEFILE=`basename ${FILE} ".bdg"`;

        SCRIPT=bdg2bw_${ID}.sh;

        echo '#!/bin/bash' >>${POUT}/${SCRIPT};
        echo '' >>${POUT}/${SCRIPT};
        echo "${PCODE}/bedtools-2.17.0/bin/slopBed -i ${FILE} -g ${PCHROM} -b 0 | ${PCODE}/UCSC_tools/bedClip stdin ${PCHROM} ${FILE}.clip" >>${POUT}/${SCRIPT};
        echo '' >>${POUT}/${SCRIPT};
        echo "${PCODE}/UCSC_tools/bedGraphToBigWig ${FILE}.clip ${PCHROM} ${POUT}/${BASEFILE}.bw" >>${POUT}/${SCRIPT};
        echo '' >>${POUT}/${SCRIPT};
        echo "rm -f ${FILE}.clip" >>${POUT}/${SCRIPT};
        echo '' >>${POUT}/${SCRIPT};

        nice -5 qsub -e ${POUT}/${ID}.err -o ${POUT}/${ID}.out -q newnodes.q ${POUT}/${SCRIPT};
        rm ${POUT}/${SCRIPT};
done

JOB_STATUS=$(qstat);
if [[ $JOB_STATUS =~ bdg2bw ]] ;
then
        COMPLETE=0;
else
        COMPLETE=1;
fi

while [ "${COMPLETE}" = "0" ]
do
        echo "not all bdg2bw jobs complete. Waiting 60 seconds.."
        sleep 60;

        JOB_STATUS=$(qstat);
        if [[ $JOB_STATUS =~ bdg2bw ]] ;
        then
                COMPLETE=0;
        else
                COMPLETE=1;
        fi
done
echo "All jobs completed, continuing to stage 4..";
echo " ";

#housekeeping
rm ${POUT}/*.bdg;
LOGS=${POUT}/d_logs;
mkdir ${LOGS};
mv "${POUT}/"*.log ${LOGS};
mv "${POUT}/"*.err ${LOGS}; 
mv "${POUT}/"*.out ${LOGS};

PBW_TREAT=${POUT}/d_treat;
PBW_LOGLR=${POUT}/d_loglr;
PBW_FE=${POUT}/d_fe;
PBW_POIS=${POUT}/d_poiss;
mkdir ${PBW_TREAT};
mkdir ${PBW_LOGLR};
mkdir ${PBW_FE};
mkdir ${PBW_POIS};
mv "${POUT}/"*treat_pileup.bw ${PBW_TREAT};
mv "${POUT}/"*logLR.bw ${PBW_LOGLR};
mv "${POUT}/"*FE.bw ${PBW_FE};
mv "${POUT}/"*ppois.bw ${PBW_POIS};

#do big wig correlate
#wigcorrelate
#/net/isi-scratch/giuseppe/tools/UCSC_tools/wigCorrelate 
#wigCorrelate - Produce a table that correlates all pairs of wigs.
#usage:
#   wigCorrelate one.wig two.wig ... n.wig
#This works on bigWig as well as wig files.
#The output is to stdout
#options:
#   -clampMax=N - values larger than this are clipped to this value

#the following needs to be done 4 times
SPACE="  ";
for i in PBW_TREAT PBW_FE PBW_POIS PBW_LOGLR
do
 	for FILE in ${!i}/*.bw;
        do
        	INPUT+=${FILE}${SPACE};
	done
	SCRIPT=wc_${i}.sh;
	echo '#!/bin/bash' >>${!i}/${SCRIPT};
	echo '' >>${!i}/${SCRIPT};
	echo "${PCODE}/UCSC_tools/wigCorrelate ${INPUT}" >>${!i}/${SCRIPT};
	echo '' >>${!i}/${SCRIPT};
	nice -5 qsub -e ${!i}/wc.err -o ${!i}/wc.out -q newnodes.q ${!i}/${SCRIPT};
	INPUT='';
done

JOB_STATUS=$(qstat);
if [[ $JOB_STATUS =~ wc ]] ;
then
        COMPLETE=0;
else
        COMPLETE=1;
fi

while [ "${COMPLETE}" = "0" ]
do
        echo "not all wc jobs complete. Waiting 60 seconds.."
        sleep 5;

        JOB_STATUS=$(qstat);
        if [[ $JOB_STATUS =~ wc ]] ;
        then
                COMPLETE=0;
        else
                COMPLETE=1;
        fi
done
echo "All jobs completed, continuing to stage 5..";
echo " ";
#--------------------------------------------------------------------------------------------------------------
#5. convert wc.out to matrix and send it to R
#--------------------------------------------------------------------------------------------------------------
#USAGE: wigcorrelate_to_R_heatmap.pl -i=<INFILE>
perl "${PSCRIPT}/"wigcorrelate_to_Rheatmap.pl -i="${PBW_TREAT}/"wc.out;
perl "${PSCRIPT}/"wigcorrelate_to_Rheatmap.pl -i="${PBW_LOGLR}/"wc.out;
perl "${PSCRIPT}/"wigcorrelate_to_Rheatmap.pl -i="${PBW_FE}/"wc.out;
perl "${PSCRIPT}/"wigcorrelate_to_Rheatmap.pl -i="${PBW_POIS}/"wc.out;

#output of the perl script is a file called wc_R.dat
#create an R script and execute it to get a pdf

SCRIPT=wc_compute.r;
for i in PBW_TREAT PBW_FE PBW_POIS PBW_LOGLR
do
	echo 'library(gplots)' >>${!i}/${SCRIPT};
	echo "vdr <- read.table(\""${!i}/"wc_R.dat\",header=TRUE, row.names=1)" >>${!i}/${SCRIPT}; 
	echo "pdf(\""${!i}/"heatmap.pdf\", height=10, width=10)" >>${!i}/${SCRIPT};
	echo 'heatmap.2(as.matrix(vdr),trace="none")' >>${!i}/${SCRIPT};
	echo 'dev.off()' >>${!i}/${SCRIPT};
	R --vanilla < ${!i}/wc_compute.r; 
done

echo "Finished.";
