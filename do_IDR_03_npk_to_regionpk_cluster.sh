#!/bin/bash

#IDR stage 3
#sort -k 8nr,8nr chipSampleRep1_VS_controlSampleRep0_peaks.encodePeak | head -n 100000 | gzip -c > chipSampleRep1_VS_controlSampleRep0.regionPeak.gz


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "<PATH> data path for the MACS2 narrowPeak files"
        exit
fi

PDATA=$1;

for FILE in ${PDATA}/*.narrowPeak;
	do ID=`basename ${FILE} ".narrowPeak"`;

        SCRIPT=encdpk_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	#echo "sort -k 8nr,8nr ${FILE} | head -n 100000 | gzip -c > ${PDATA}/${ID}.regionPeak.gz"
	echo "sort -k 8nr,8nr ${FILE} | head -n 100000 > ${PDATA}/${ID}.regionPeak"  >>${PDATA}/${SCRIPT};
        
	nice -5 qsub -e ${PDATA}/encdpk_${ID}.err  -o ${PDATA}/encdpk_${ID}.out  -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done

