#!/bin/bash

#wrapper to run do_pscanchip_out_intersect_vdrbvs.pl on the cluster

#Input: a directory full of .ris files, output of PscanCHIP. each for a different TF whose PWM intervals around your sequences you've chosen to download
#This script will output a SUBSET of the input intervals, including intervals WHICH INTERSECT a VDR-BV list provided in input
#You can also threshold by pscanchip score.


if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <VDR-BV_VCF> <PSCANCHIP_THRS>"
        echo "<PATH> - absolute path for the pscanchip.ris files, one per TF"
	echo "<VDR-BV_VCF> vcf containing the VDR-BVs or VDR-rBVs, full path";
	echo "<PSCANCHIP_THRS> [0.0-1.0]  NA if no threshold";
        exit
fi

PDATA=$1;
VCFDATA=$2;
THRS=$3
PCODE="/net/isi-backup/giuseppe/scripts/P_Various/do_pscanchip_out_intersect_vdrbvs.pl";

for FILE in ${PDATA}/Pscanchip*.ris;
        do
#	ID=`echo ${FILE} | grep -Po "BVs_(.*)_sites"`;
        SEED=`echo ${FILE} | grep -Po "BVs_(.*)_sites"`;
        TF=`echo ${SEED} | cut -d _ -f 2`;
        ID=`echo ${SEED} | cut -d _ -f 3`;

	BASENAME=`basename ${FILE} ".ris"`;


	if [ $THRS == 'NA' ]; then 
		SCRIPT=script_${ID}.sh;
	else
		SCRIPT=script_m${THRS}_${ID}.sh;
	fi

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        if [ $THRS == 'NA' ]; then
		echo "perl ${PCODE} -v=${VCFDATA} -p=${FILE}" >>${PDATA}/${SCRIPT};
	else
		echo "perl ${PCODE} -v=${VCFDATA} -p=${FILE} -m=${THRS} " >>${PDATA}/${SCRIPT};
	fi

	if [ $THRS == 'NA' ]; then
		nice -5 qsub -e ${PDATA}/script_${ID}.err -o ${PDATA}/Pscanchip_allsites_BVs_${TF}_${ID}_sites.out -q medium_jobs.q ${PDATA}/${SCRIPT};
	else
		nice -5 qsub -e ${PDATA}/script_m${THRS}_${ID}.err -o ${PDATA}/Pscanchip_sites_m${THRS}_BVs_${TF}_${ID}_sites.out -q medium_jobs.q ${PDATA}/${SCRIPT};
	fi
        rm ${PDATA}/${SCRIPT};  
done
