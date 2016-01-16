#!/bin/bash

#Script to create two "mutant" peaklist bed files for MAnorm.
#The idea is the following
#Manorm builds its model from the COMMON peaks from peakset1 and peakset2
#but we want a better way to decide what is a common peak
#Separately,I had obtained beds of common peaks using the IDR method

#in this script, I will create 2 peaksets as follows
#both peaksets: IDR peaks
#peakset 1: given a bed of called peaks (say GEM) add the difference between GEM and IDR peaks
#peakset 2: given a bed of called peaks (say GEM) add the difference between GEM and IDR peaks

#INPUT
#1 bed IDR overlapping peaks
#2.1 bed GEM(MACS) peaks sample 1
#2.2 bed GEM(MACS) peaks sample 2



if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the bed files"
        exit
fi

PDATA=$1;

#awk -v t=$2 '$11 <= t {print $0}' $1

for FILE in ${PDATA}/*.bed;
        do
        BASEFILE=`basename ${FILE} ".bed"`;
        SCRIPT=MANORM_${BASEFILE}.sh;

	ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        #change required options here:
        echo "awk '{print \$1, \"\t\", \$2, \"\t\", \$3}' ${FILE}" >>${PDATA}/${SCRIPT};

        nice -5 qsub -o ${PDATA}/${ID}_MAnorm_peaks.bed -q newnodes.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT};
done
#for FILE in ${PDATA}/*.bed;
#	do 
#        BASEFILE=`basename ${FILE} ".bed"`;  
#	awk '{print $1, "\t", $2, "\t", $3, "\t", $6}' ${FILE} > ${BASEFILE}_manorm.bed
#done
