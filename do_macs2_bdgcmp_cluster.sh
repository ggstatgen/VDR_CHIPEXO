#!/bin/bash

#scripts to do the following:
#macs2 bdgcmp -t GM19249_treat_pileup.bdg -c GM19249_control_lambda.bdg -o GM19249_ppois.bdg                     >& bdgcmp_GM..._ppois.out    &
#macs2 bdgcmp -t GM19249_treat_pileup.bdg -c GM19249_control_lambda.bdg -o GM19249_ko_FE.bdg -m FE               >& bdgcmp_GM_FE.out    &
#macs2 bdgcmp -t GM19249_treat_pileup.bdg -c GM19249_control_lambda.bdg -o GM19249_logLR.bdg -m logLR -p 0.00001 >& bdgcmp_GM_logLR.out &



if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "You must specify the data path for the bdg files (e.g. /home/me/files)"
        exit
fi

PDATA=$1;
#location of macs2 executable, change if needed
PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";

for FILE in ${PDATA}/*.bdg;
	#do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
	do ID=`echo ${FILE} | egrep -o "[0-9]D3*"`;
	echo $ID;
	#do ID=`basename ${FILE} ".bdg"`;

	C=`echo ${FILE} | egrep -o "control"`;
        if [ "${C}" = 'control' ];
        then
                continue;
        fi

        SCRIPT=macs2_bdgcmp_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
        echo "${PCODE}/macs2 bdgcmp -t ${PDATA}/${ID}_treat_pileup.bdg -c ${PDATA}/${ID}_control_lambda.bdg -o ${PDATA}/${ID}_ppois.bdg" >>${PDATA}/${SCRIPT};
        echo "${PCODE}/macs2 bdgcmp -t ${PDATA}/${ID}_treat_pileup.bdg -c ${PDATA}/${ID}_control_lambda.bdg -o ${PDATA}/${ID}_FE.bdg -m FE" >>${PDATA}/${SCRIPT};
        echo "${PCODE}/macs2 bdgcmp -t ${PDATA}/${ID}_treat_pileup.bdg -c ${PDATA}/${ID}_control_lambda.bdg -o ${PDATA}/${ID}_logLR.bdg -m logLR -p 0.00001" >>${PDATA}/${SCRIPT};
        nice -5 qsub -e ${PDATA}/macs2_${ID}.err -o ${PDATA}/macs2_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        #rm ${PDATA}/${SCRIPT};  
done
