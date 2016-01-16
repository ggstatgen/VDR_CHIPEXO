#!/bin/bash


#run wigcorrelated with on all the bigwig in a directory

#wigcorrelate
#/net/isi-scratch/giuseppe/tools/UCSC_tools/wigCorrelate 

#wigCorrelate - Produce a table that correlates all pairs of wigs.
#usage:
#   wigCorrelate one.wig two.wig ... n.wig
#This works on bigWig as well as wig files.
#The output is to stdout
#options:
#   -clampMax=N - values larger than this are clipped to this value

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input .bw files"
        exit
fi

PDATA=$1;
PCODE='/net/isi-scratch/giuseppe/tools/UCSC_tools';

for FILE in ${PDATA}/*.bw;
        do
        #OPT="-I ";
        SPACE="  ";
        INPUT+=${FILE}${SPACE};
        #need to build list of files
        #-I file1 -I file2 -I file3
done

nice nohup ${PCODE}/wigCorrelate ${INPUT} >& ${PDATA}/wigcorrelate.out  &

#for FILE in ${PDATA}/*.wg;
#   do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
#   ${PCODE}macs14 -t ${FILE} -n ${PDATA}${ID} -g hs --keep-dup=all -w -S;
#done
