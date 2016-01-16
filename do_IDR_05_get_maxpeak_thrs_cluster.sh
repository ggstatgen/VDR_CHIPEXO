#!/bin/bash

#Run this on the output of the idr


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "<PATH> data path for the -overlapped-peaks.txt directory"
	echo "NOTE IDR - 0.1"
        exit
fi

PDATA=$1;

for FILE in ${PDATA}/*-overlapped-peaks.txt;
	do ID=`basename ${FILE} "-overlapped-peaks.txt"`;

        SCRIPT=idr_pkthrs_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "awk '\$11 <= 0.1 {print \$0}' ${FILE} | wc -l"  >> ${PDATA}/${SCRIPT};
        
	nice -5 qsub -e ${PDATA}/idr_pkthrs_${ID}.err  -o ${PDATA}/idr_pkthrs_${ID}.out  -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done

#cat ${PDATA}/idr_pkthrs_*.out > ${PDATA}/IDR_peak_thresholds.out
