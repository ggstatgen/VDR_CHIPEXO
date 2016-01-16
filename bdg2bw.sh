#!/bin/bash

CODE="/net/isi-scratch/giuseppe/tools";

# check commands: slopBed, bedGraphToBigWig and bedClip

which slopBed &>/dev/null || which bedtools &>/dev/null || { echo "slopBed/bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }

# end of checking

if [ $# -lt 2 ];then
    echo "Need 2 parameters! <bedgraph> <chrom info>"
    exit
fi

F=$1
G=$2

${CODE}/bedtools-2.17.0/bin/slopBed -i ${F} -g ${G} -b 0 | ${CODE}/UCSC_tools/bedClip stdin ${G} ${F}.clip

${CODE}/UCSC_tools/bedGraphToBigWig ${F}.clip ${G} ${F/bdg/bw}

rm -f ${F}.clip
