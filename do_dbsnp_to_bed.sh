#!/bin/bash

#found this here
#https://gist.github.com/arq5x/859487#file-dbsnp-to-bed-sh

export GENOME=hg19
export SNPBUILD=138

curl -s http://hgdownload.cse.ucsc.edu/goldenPath/$GENOME/database/snp$SNPBUILD.txt.gz | \
zcat | \
cut -f 2,3,4,5,6,7,10,16 > dbsnp.$SNPBUILD.$GENOME.bed
 
#head dbsnp.$SNPBUILD.$GENOME.bed
#chr1 10433 10433 rs56289060 0 + -/C near-gene-5
#chr1 10491 10492 rs55998931 0 + C/T near-gene-5
#chr1 10518 10519 rs62636508 0 + C/G near-gene-5
#chr1 10582 10583 rs58108140 0 + A/G near-gene-5
#chr1 10827 10828 rs10218492 0 + A/G near-gene-5
#chr1 10903 10904 rs10218493 0 + A/G near-gene-5
#chr1 10926 10927 rs10218527 0 + A/G near-gene-5
#chr1 10937 10938 rs28853987 0 + A/G near-gene-5
#chr1 11001 11002 rs79537094 0 + A/C near-gene-5
#chr1 11013 11014 rs28484712 0 + A/G near-gene-5
