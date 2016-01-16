#!/bin/bash

#scripto to run macs2 on tagalign files for the IDR pipeline
# https://sites.google.com/site/anshulkundaje/projects/idr

#hg19 mappability sizes
#from https://groups.google.com/forum/?fromgroups=#!topic/macs-announcement/9DeYzN3vFWA
# 40	2540757438

#mm9
#51 2271970228

#macs2 callpeak -t chipSampleRep1.tagAlign.gz -c controlSampleRep0.tagAlign.gz -f BED -n chipSampleRep1_VS_controlSampleRep0 -g hs -p 1e-3 --to-large

if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <PATH> <EXT_SIZE> <DUPLICATES> <PVAL>"
        echo "<PATH> data path for the tagAlign.gz files"
        echo "<EXT_SIZE> extension size in bp"
        echo "<DUPLICATES> one of [1|auto|all]"
	echo "<PVAL> p value threshold, eg 1e-3"
        exit
fi

PDATA=$1;
EXT_SIZE=$2;
DUPS=$3;
PVAL=$4;
#location of macs2 executable, change if needed
PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";

for FILE in ${PDATA}/*.tagAlign.gz;
	do ID=`basename ${FILE} ".tagAlign.gz"`;

        SCRIPT=MACS2_d${EXT_SIZE}_dup${DUPS}_p${PVAL}_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "source activate" >> ${PDATA}/${SCRIPT};

	echo "${PCODE}/macs2 callpeak -t ${FILE} -f BED -n  ${PDATA}/MACS2_d${EXT_SIZE}_dup${DUPS}_p${PVAL}_${ID} --keep-dup=${DUPS}  -g 2540757438 -p ${PVAL} --nomodel --extsize ${EXT_SIZE}" >>${PDATA}/${SCRIPT};


        #echo "${PCODE}/macs2 callpeak -t ${FILE} -f BED -n  ${PDATA}/MACS2_d${EXT_SIZE}_dup${DUPS}_p${PVAL}_${ID} --keep-dup=${DUPS}  -g 2540757438 -p ${PVAL} --nomodel --extsize ${EXT_SIZE}  --to-large" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/MACS2_d${EXT_SIZE}_dup${DUPS}_p${PVAL}_${ID}.err -o  ${PDATA}/MACS2_d${EXT_SIZE}_dup${DUPS}_p${PVAL}_${ID}.out  -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
