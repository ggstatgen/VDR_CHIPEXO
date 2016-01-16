#!/bin/bash


#get best n peaks by p-value (column 8) in an extended bed file (narrowpeak etc) from ENCODE specs

if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <NUMBER> <EXT>"
        echo "<PATH> data path for the interval files"
	echo "<NUMBER> top peaks to keep"
	echo "<EXT_SIZE> file extension"
        exit
fi

PDATA=$1;
NUMBER=$2;
EXT=$3;

for FILE in ${PDATA}/*.${EXT};
	do ID=`basename ${FILE} ".${EXT}"`;

        SCRIPT=npk_get_best_${NUMBER}_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "cat ${FILE} | sort -k8nr,8nr | head -n ${NUMBER} | sort -k1,1V -k2,2g" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/npk_get_best_${NUMBER}_${ID}.err -o ${PDATA}/best_${NUMBER}_${ID}.peaks  -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
