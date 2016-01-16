#!/bin/bash

#script to "distil" manorm compatible-BED from the bed peaks created by GPS and MACS
#input: chr1    10239   10279   other other ....
#output: chr1     10239   10279
#command: awk '{print $1, "\t", $2, "\t", $3}'  GM06986_1_GEM_events.bed > GM06986_peaks.bed

#GPS has a header, you might need to try to remove it automatically with something like head -n1


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
