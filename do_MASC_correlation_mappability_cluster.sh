#!/bin/bash

#MaSC: Mappability-Sensitive Cross-Correlation for Estimating Fragment Length

#Paper
#http://bioinformatics.oxfordjournals.org/content/29/4/444

#explanation
#https://www.biostars.org/p/63720/

#software download
#http://www.perkinslab.ca/Software.html

#This should do the same as Anshul's software phantompeak tool, but improve things because it accounts for mappability

#synopsis
#Usage:
#    MaSC.pl --help --verbose --mappability_path=/mappability/dir/
#    --chrom_length_file=/dir/chrom_lens.txt --input_bed=/dir/input.bed
#    --prefix=myprefix --smooth_win_size=15 --min_shift=0 --max_shift=400


if [ ! $# == 5 ]; then
        echo "Usage: `basename $0` <PATH> <MAP_PATH> <MINSHIFT> <MAXSHIFT> <SMOOTH_WIN>"
        echo "<PATH> for the BED ALIGNMENT files"
	echo "<MAP_PATH> path of the directory containing MASC mappability wigs"
	echo "<MINSHIFT> eg 0"
	echo "<MAXSHIFT> eg 400"
	echo "<SMOOTH_WIN> eg 15"
	echo "NOTE: edit script to change paths to mappability files. Currently hg19, 40bp"
        exit
fi

PPERL="/net/isi-cgat/ifs/apps/apps/perl-5.16.1/bin";
PDATA=$1;
PMAPPABILITY=$2;
MINSHIFT=$3;
MAXSHIFT=$4;
SMOOTH=$5;
PCODE="/net/isi-scratch/giuseppe/tools/MASC/MaSC.pl";
#PMAPPABILITY="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/MAPPABILITY/MASC/"
#PCHROMINFO="/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes"
PCHROMINFO="/net/isi-scratch/giuseppe/indexes/chrominfo/mm10.chrom.sizes"
mm10.chrom.sizes


for FILE in ${PDATA}/*.bed;
	do ID=`basename ${FILE} ".bed"`;

        SCRIPT=MaSC_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "${PPERL}/perl ${PCODE} --verbose --mappability_path=${PMAPPABILITY} --chrom_length_file=${PCHROMINFO} --input_bed=${FILE} --prefix=${PDATA}/MaSC_${ID}  --smooth_win_size=${SMOOTH} --min_shift=${MINSHIFT} --max_shift=${MAXSHIFT}" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/MaSC_${ID}.err -o ${PDATA}/MaSC_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
