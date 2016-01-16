#!/bin/bash

#Giuseppe 15/4/2013
#Taken from idr_count_peaks.sh
#this returns all the peaks which pass the threshold
#you need the output of this file to do the FRIP and marginal enrichment calculations

#taken from:
#Script to count peaks after IDR analysis
#from https://sites.google.com/site/anshulkundaje/projects/idr#TOC-Summary
#THRESHOLD
#Anshul suggests a threshold of 0.1 for preudoreplicates with shallow depth:
#For self-consistency analysis of datasets with shallow sequencing depths, you can use an IDR threshold as relaxed as 0.1 if you start with < 100K pre-IDR peaks.

if [ $# -lt 2 ];then
    echo "USAGE: idr_count_peaks.sh <FILENAME> <THRESHOLD>"
    echo "FILENAME: idr overlap file"
    echo "THRESHOLD: idr threshold (e.g. 0.01)"
    exit
fi

awk -v t=$2 '$11 <= t {print $0}' $1


