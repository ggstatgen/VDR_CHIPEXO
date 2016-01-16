#!/usr/bin/env python

import os

os.system('bamToBed -i something.sorted.bam > something.sorted.bed')
os.system('genetrack something.sorted.bed -b > something.sorted.peaks.gff')

slopBed -i forward.bedgraph -g sacCer3_size.txt -l 0 -r 200 > forward.slop.bedgraph

intersectBed -a forward.slop.bedgraph -b reverse.bedgraph -wo -bed > intersected.txt

awk -f tobed.awk intersect.txt > fullpeaks.bed




#for example
bedtools slop -i FLFN_F1804_thymus_forward.bedgraph -g /net/isi-scratch/giuseppe/indexes/chrominfo/mm10.chrom.sizes -l 0 -r 100 > FLFN_F1804_thymus_forward_slop.bedgraph
bedtools intersect -a FLFN_F1804_thymus_forward_slop.bedgraph -b FLFN_F1804_thymus_reverse.bedgraph -wo -bed > intersected.txt
