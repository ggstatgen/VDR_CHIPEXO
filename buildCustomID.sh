#!/bin/bash

#update "nomono" version
for chr in $(seq 1 22); do
	zcat ALL.chr$chr\.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz | awk -v c=$chr '{ if ($1 == ".") { lr=length($3); la=length($4); if (lr == 1 && la == 1) $1="chr"c":"$2":S"; else if (lr > la) $1="chr"c":"$2":D"; else $1="chr"c":"$2":I";} print $0 }' | gzip -c > ALL.chr$chr\.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.withCustomID.legend.gz
done

#update "nosing" version
#for chr in $(seq 1 22); do
#	zcat ALL.chr$chr\.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz | awk -v c=$chr '{ if ($1 == ".") { lr=length($3); la=length($4); if (lr == 1 && la == 1) $1="chr"c":"$2":S"; else if (lr > la) $1="chr"c":"$2":D"; else $1="chr"c":"$2":I";} print $0 }' | gzip -c > ALL.chr$chr\.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.withCustomID.legend.gz
#done


