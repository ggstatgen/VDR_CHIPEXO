#!/bin/bash

#script to run macs2 with custom options on the cluster

#hg19 mappability sizes
#from https://groups.google.com/forum/?fromgroups=#!topic/macs-announcement/9DeYzN3vFWA
# 40	2540757438

#mm9
#51 2271970228

# 19/12/2012 
# VDR data: I could probably use genetrack's estimate for the peak pair distance as an input to macs2!!
# genetrack has no statistical modelling of the background, however in the document coming with the data (pag 5)
# it says that c-w distance (bp) is 26 : "base pair distance between paired forward strand peak and reverse strand peak"

if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <PATH> <EXT_SIZE> <DUPLICATES> <QVAL>"
        echo "<PATH> data path for the bam files"
	echo "<EXT_SIZE> extension size in bp"
	echo "<DUPLICATES> one of [1|auto|all]"
	echo "<QVAL> q value threshold, eg 0.01"
        exit
fi

PDATA=$1;
EXT_SIZE=$2;
DUPS=$3;
QVAL=$4;
#location of macs2 executable, change if needed
PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";

for FILE in ${PDATA}/*.bam;
	do ID=`basename ${FILE} ".bam"`;

        SCRIPT=MACS2_d${EXT_SIZE}_dup${DUPS}_q${QVAL}_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "source activate" >> ${PDATA}/${SCRIPT};

	#change required options here:
	echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/MACS2_${ID}_d${EXT_SIZE}_dups_${DUPS}_q_${QVAL} -g mm --keep-dup=${DUPS} --nomodel --extsize ${EXT_SIZE} -q 0.01" >>${PDATA}/${SCRIPT};
	#echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/MACS2_${ID}_d${EXT_SIZE}_dups_${DUPS}_q_${QVAL} -g mm" >>${PDATA}/${SCRIPT};

	#echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/MACS2_${ID} -g 2540757438 --keep-dup=auto" >>${PDATA}/${SCRIPT};

        #echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/MACS2_d${EXT_SIZE}_dup${DUPS}_q${QVAL}_${ID} -g 2540757438 --keep-dup=${DUPS} -q ${QVAL}  -B --SPMR  --nomodel --extsize ${EXT_SIZE}" >>${PDATA}/${SCRIPT};
	#echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/MACS2_${ID} -g 2540757438 -B --SPMR --keep-dup=all -q 0.01 --nomodel --extsize 13" >>${PDATA}/${SCRIPT};
        #echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/${ID} -g 2540757438 -B --trackline --SPMR --keep-dup=all -s 40 -p 1e-5 --nomodel --extsize 26" >>${PDATA}/${SCRIPT};
        
        #extra loose thresholds FOR generating GATK CONSENSUS INTERVALS ONLY
        #echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/MACS2_${ID} -g 2540757438 --keep-dup=1  -B --SPMR -q 0.01 --nomodel --extsize 26  --call-summits" >> ${PDATA}/${SCRIPT};
	#echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/MACS2_${ID} -g 2540757438 --keep-dup=auto --nomodel --extsize 26 -q 0.000000001" >> ${PDATA}/${SCRIPT};
	#echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/${ID} -g 2540757438 --keep-dup=auto -s 40 --nomodel --extsize 26" >>${PDATA}/${SCRIPT};        
	#echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/MACS2_${ID} -g 2540757438 --keep-dup=auto -s 40 -q 0.01 --broad --nomodel --extsize 26" >>${PDATA}/${SCRIPT};
	#echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/MACS2_${ID} -g mm --keep-dup=all --nomodel --extsize 39" >> ${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/MACS2_d${EXT_SIZE}_dup${DUPS}_q${QVAL}_${ID}.err -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
