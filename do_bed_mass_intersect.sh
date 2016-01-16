#!/bin/bash

#Intersect all the beds in a directory with another bed
#I need it to intersect all the peaks from my chipexo samples with the peaks in the chipseq consensus

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <BED_FILE"
        echo "<PATH> directory for the .bed files"
	echo "<BED_FILE> full path to file you want to intersect with"
	echo "NOTE: MAKE SURE the genome reference is the same!"
        exit
fi

PDATA=$1;
PBED=$2;
PCODE="/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin/bedtools";
#GENOME="/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes";
#GENOME="/net/isi-scratch/giuseppe/indexes/chrominfo/mm10.chrom.sizes";

for FILE in ${PDATA}/*.bed;
	do ID=`basename ${FILE} ".bed"`;

        SCRIPT=mass_intr_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	
	echo "source activate" >> ${PDATA}/${SCRIPT};
	echo "${PCODE} intersect -a ${FILE}  -b ${PBED}" >> ${PDATA}/${SCRIPT};
	#echo "${PCODE} intersect -b ${PSLOP} -i ${FILE} -g ${GENOME} | sort -k1,1V -k2,2g | ${PCODE} merge -i stdin" >> ${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/mass_intr_${ID}.err -o ${PDATA}/${ID}_intersect_CHIPSEQ.bed -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
