#!/bin/bash

#script to "distil" manorm compatible-BED from the bed alignment extracted by bedtools bamtobed
#input: chr1    10239   10279   DGM97JN1_120615_0211_AC0R2FACXX:3:2106:17704:25973#0    255     -
#output: chr1     10239   10279   -
#command: awk '{print $1, "\t", $2, "\t", $3, "\t", $6}' in.bed > out.bed

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
        SCRIPT=${BASEFILE}.sh;

	ID=`echo ${FILE} | egrep -o "GM[0-9]*"`; 
	#dataset dependent

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        #change required options here:
        echo "awk '{print \$1, \"\t\", \$2, \"\t\", \$3, \"\t\", \$6}' ${FILE}" >>${PDATA}/${SCRIPT};

        nice -5 qsub -o ${PDATA}/${ID}_MAnorm_reads.bed -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done



#for FILE in ${PDATA}/*.bed;
#	do 
#        BASEFILE=`basename ${FILE} ".bed"`;  
#	awk '{print $1, "\t", $2, "\t", $3, "\t", $6}' ${FILE} > ${BASEFILE}_manorm.bed
#done
