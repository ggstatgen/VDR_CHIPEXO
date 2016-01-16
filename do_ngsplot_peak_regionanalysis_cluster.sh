#!/bin/bash

#This runs regionanalysis
#https://github.com/shenlab-sinai
#to see which intervals in an input bed file overlap genomic features (annotate peaks with genomic features)
#It is probably useful as a preprocessing for ngs_analysis. You can select only genes which overlap peaks, ano not all genes.
#You can then do multiple plots to see how different sets of peak/gene behave compared to all genes

#given the output of this, I generate genelists for ngsplot using
#grep -P "romoter" vdr_o10_peaks.bed.annotated | cut -f 7 | sort | uniq > vdr_o10_peaks_promoters.list
#grep -P "body" vdr_o10_peaks.bed.annotated | cut -f 7 | sort | uniq > vdr_o10_peaks_promoters.list

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_BED> <DB>"
	echo "<PATH_BED> path for the bed files with peaks to annotate"
	echo "<DB> database to use [ensembl|refseq]"
        exit
fi

PDATA=$1;
PDB=$2;
PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";
GENOME="hg19";
		
for FILE in ${PDATA}/*.bed;
	do ID=`basename ${FILE} ".bed"`;

        SCRIPT=region_a_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate'  >>${PDATA}/${SCRIPT};
	echo "${PCODE}/region_analysis.py -i ${FILE} -d ${PDB} -g ${GENOME}" >> ${PDATA}/${SCRIPT};

        nice -5 qsub -e ${POUT}/region_a_${ID}.err -o ${POUT}/region_a_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
