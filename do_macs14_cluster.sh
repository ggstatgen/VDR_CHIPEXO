#!/bin/bash

#script to run macs14 with custom options on the cluster

#hg19 mappability sizes
#from https://groups.google.com/forum/?fromgroups=#!topic/macs-announcement/9DeYzN3vFWA
# 40    2540757438


#mouse 
# 51  2271970228


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
#location of fastqc executable, change if needed
PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";

for FILE in ${PDATA}/*.bam;
	#do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
	do ID=`basename ${FILE} ".bam"`;
        #do ID=`echo ${FILE} | egrep -o "38380_[0-9]*"`;        

        SCRIPT=MACS14_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
        #echo "${PCODE}/macs14 -t ${FILE} -n ${PDATA}/${ID} -s 40 --nomodel --shiftsize=13 --keep-dup=all -p 1e-8 -g 2540757438" >>${PDATA}/${SCRIPT};
        echo "${PCODE}/macs14 -t ${FILE} -n ${PDATA}/MACS14_${ID} -s 40 --nomodel --shiftsize=13 --keep-dup=1 -g 2540757438 -w -S" >>${PDATA}/${SCRIPT};
        
        nice -5 qsub -e ${PDATA}/MACS14_${ID}.err -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done


#for FILE in ${PDATA}*.bam;
#   do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
#   ${PCODE}/macs14 -t ${FILE} -n ${PDATA}${ID} -g hs --keep-dup=all -w -S;
#done
