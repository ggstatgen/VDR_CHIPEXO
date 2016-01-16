#!/bin/bash

#Tool:    bedtools slop (aka slopBed)
#Version: v2.17.0
#Summary: Add requested base pairs of "slop" to each feature.
#
#Usage:   bedtools slop [OPTIONS] -i <bed/gff/vcf> -g <genome> [-b <int> or (-l and -r)]
#
#Options: 
#	-b	Increase the BED/GFF/VCF entry -b base pairs in each direction.
#		- (Integer) or (Float, e.g. 0.1) if used with -pct.
#
#	-l	The number of base pairs to subtract from the start coordinate.
#		- (Integer) or (Float, e.g. 0.1) if used with -pct.
#
#	-r	The number of base pairs to add to the end coordinate.
#		- (Integer) or (Float, e.g. 0.1) if used with -pct.
#
#	-s	Define -l and -r based on strand.
#		E.g. if used, -l 500 for a negative-stranded feature, 
#		it will add 500 bp downstream.  Default = false.
#
#	-pct	Define -l and -r as a fraction of the feature's length.
#		E.g. if used on a 1000bp feature, -l 0.50, 
#		will add 500 bp "upstream".  Default = false.
#
#	-header	Print the header from the input file prior to results.


if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <SLOP> <EXT>"
        echo "<PATH> directory for the interval files"
	echo "<SLOP> amount to extend left and right"
	echo "<EXT> type of interval (eg gff, bed)"
        exit
fi

PDATA=$1;
PSLOP=$2;
EXT=$3;
PCODE="/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin/bedtools";
GENOME="/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes";
#GENOME="/net/isi-scratch/giuseppe/indexes/chrominfo/mm10.chrom.sizes";

for FILE in ${PDATA}/*.${EXT};
	do ID=`basename ${FILE} ".${EXT}"`;

        SCRIPT=slopmerge_${PSLOP}_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	
	echo "source activate" >> ${PDATA}/${SCRIPT};
	echo "${PCODE} slop -b ${PSLOP} -i ${FILE} -g ${GENOME} | sort -k1,1V -k2,2g | ${PCODE} merge -i stdin" >> ${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/slopmerge_${PSLOP}_${ID}.err -o ${PDATA}/${ID}_slop${PSLOP}merge.${EXT} -q newnodes.q ${PDATA}/${SCRIPT};
        #rm ${PDATA}/${SCRIPT};  
done
