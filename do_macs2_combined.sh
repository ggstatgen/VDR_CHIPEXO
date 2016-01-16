#!/bin/bash

#script to run macs2 on all bam files in a directory
#follows the ideas on this page
#https://github.com/taoliu/MACS/wiki/Build-Signal-Track
#Ideally, you have looked at the correlation between some wigs as specified
#Then, choose the files which correlate well and run this script on them

#hg19 mappability sizes
#from https://groups.google.com/forum/?fromgroups=#!topic/macs-announcement/9DeYzN3vFWA
# 40	2540757438

if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <MACSDUP> <MACS_Q>"
        echo "PATH - path for the bam files (e.g. /home/me/files)"
        echo "MACSDUP - one of [all|auto|1]"
        echo "MACS_Q - eg 0.001 (3) or 0.00001 (5)"
        exit
fi

PDATA=$1;
MACSDUP=$2; #MACSDUP=all;
MACSQ=$3; #MACSQ=0.00001;
MAPPABILITY="2540757438";
#location of macs2 executable, change if needed
PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";

#create output dir
POUT=${PDATA}/d_MACS2_comb_${MACSDUP}_q${MACSQ};
mkdir ${POUT};

SPACE=" ";
for FILE in ${PDATA}/*.bam;
        do
        INPUT+=${FILE}${SPACE};
done

SCRIPT=MACS2_combo.sh;
echo '#!/bin/bash' >>${POUT}/${SCRIPT};
echo '' >>${POUT}/${SCRIPT};
echo "${PCODE}/macs2 callpeak -t ${INPUT} -n ${POUT}/VDR_macs2_combo -g ${MAPPABILITY} --keep-dup=${MACSDUP} -B --SPMR -q ${MACSQ} --nomodel --extsize 26" >> ${POUT}/${SCRIPT};
#echo "${PCODE}/macs2 callpeak -t ${INPUT} -n ${POUT}/VDR_macs2_combo -g ${MAPPABILITY} --keep-dup=${MACSDUP} -q ${MACSQ} --nomodel --extsize 26" >> ${POUT}/${SCRIPT};
#echo "${PCODE}/macs2 callpeak -t ${INPUT} -n ${POUT}/VDR_macs2_combo -g ${MAPPABILITY} --keep-dup=${MACSDUP} -q ${MACSQ}" >> ${POUT}/${SCRIPT};

nice -5 qsub -e ${POUT}/macs2_combo.err -o ${POUT}/macs2_combo.out -q newnodes.q ${POUT}/${SCRIPT};
rm ${POUT}/${SCRIPT};
