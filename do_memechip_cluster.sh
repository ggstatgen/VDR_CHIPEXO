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

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "You must specify the data path for the bam files (e.g. /home/me/files)"
        exit
fi

PDATA=$1;
#location of macs2 executable, change if needed
PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";

for FILE in ${PDATA}/*.bam;
	#do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
	do ID=`basename ${FILE} ".bam"`;
	#do ID=`echo ${FILE} | egrep -o "38380_[0-9]*"`;        

        SCRIPT=MACS2_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
        #echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/${ID} -g 2540757438 --keep-dup=all -q 0.0000000001 --nomodel --extsize 26" >>${PDATA}/${SCRIPT};
	#echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/${ID} -g 2540757438 --keep-dup=1 -q 0.01 --nomodel --extsize 26" >>${PDATA}/${SCRIPT};
        #echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/${ID} -g 2540757438 -B --trackline --SPMR --keep-dup=all -s 40 -p 1e-5 --nomodel --extsize 26" >>${PDATA}/${SCRIPT};
        
        #extra loose thresholds FOR generating GATK CONSENSUS INTERVALS ONLY
        echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/MACS2_${ID} -g 2540757438 --keep-dup=1  -B --SPMR -q 0.0000001 --nomodel --extsize 26" >> ${PDATA}/${SCRIPT};
	#echo "${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/${ID} -g 2540757438 --keep-dup=all -B -s 40 --nomodel -q 0.000000001 --extsize 26" >>${PDATA}/${SCRIPT};        

        nice -5 qsub -e ${PDATA}/MACS2_${ID}.err -o ${PDATA}/MACS2_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
