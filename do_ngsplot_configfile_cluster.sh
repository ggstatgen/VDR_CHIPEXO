#!/bin/bash

#This assumes you will specify a "config*.txt" configuration file as follows
#see https://code.google.com/p/ngsplot/wiki/HowToUseConfiguration

#1.bam<tab>-1<tab>"title"<>fragment_length
#2.bam<tab>-1<tab>"title"<>fragment_length

#This is good if you want to plot 2 bams in one file, eg chipexo and chipseq


#This runs the ngsplot program
#https://code.google.com/p/ngsplot/
#to obtain various plots of the pileup at the TSS, TTS, full gene, plus also heatmaps

#Usage: ngs.plot.r -G genome -R region -C [cov|config]file
#  -O name [Options]
#
### Mandatory parameters:
#  -G   Genome name. Use ngsplotdb.py list to show available genomes.
#  -R   Genomic regions to plot: tss, tes, genebody, exon, cgi, enhancer, dhs or bed
#  -C   Indexed bam file or a configuration file for multiplot
#  -O   Name for output: multiple files will be generated
### Optional parameters related to configuration file:
#  -E   Gene list to subset regions OR bed file for custom region
#  -T   Image title
### Important optional parameters:
#  -F   Further information provided to select database table or plottype:
# This is a string of description separated by comma.
# E.g. protein_coding,K562,rnaseq(order of descriptors does not matter)
#      means coding genes in K562 cell line drawn in rnaseq mode.
#  -D   Gene database: ensembl(default), refseq
#  -I   Shall interval be larger than flanking in plot?(0 or 1, default=automatic)
#  -L   Flanking region size
#  -N   Flanking region factor(will override flanking size)
#  -S   Randomly sample the regions for plot, must be:(0, 1]
#  -P   #CPUs to use. Set 0(default) for auto detection
### Misc. parameters:
#  -GO  Gene order algorithm used in heatmaps: total(default), hc, max,
# prod, diff, pca and none(according to gene list supplied)
#  -AL  Algorithm used to normalize coverage vectors: spline(default), bin
#  -CS  Chunk size for loading genes in batch(default=100)
#  -FL  Fragment length used to calculate physical coverage(default=150)
#  -MQ  Mapping quality cutoff to filter reads(default=20)
#  -SE  Shall standard errors be plotted?(0 or 1)
#  -RB  The fraction of extreme values to be trimmed on both ends
# default=0, 0.05 means 5% of extreme values will be trimmed
#  -RZ  Remove all zero profiles in heatmaps(default=1). Set 0 to keep them.
#  -SC  Color scale used to map values to colors in a heatmap.
# local(default): base on each individual heatmap
# region: base on all heatmaps belong to the same region
# global: base on all heatmaps together
# min_val,max_val: custom scale using a pair of numerics
#  -FC  Flooding fraction:[0, 1), default=0.02
#  -FI  Forbid image output if set to 1(default=0)
#  -MW  Moving window width to smooth avg. profiles, must be integer
# 1=no(default); 3=slightly; 5=somewhat; 9=quite; 13=super.
#  -H   Opacity of shaded area, suggested value:[0, 0.5]
# default=0, i.e. no shading, just curves


if [ ! $# == 6 ]; then
echo "Usage: `basename $0` <D_PATH> <CONFIG_FILE> <DB> <SHADING> <PROCS> <TITLE>"
	echo "<D_PATH> where the config file is and where you want the output"
	echo "<CONFIG_FILE> configuration file name in the path above"
	echo "<DB> database to use [ensembl|refseq]"
	echo "<SHADING> opacity under curve [0,..., 0.5]"
	echo "<PROCS> processors to use, 0 means all"
	echo "<TITLE> plot basename title (eg bam1 vs bam2)"
exit
fi

PDATA=$1;
CONFIG_FILE=$2;
PDB=$3;
SHADING=$4;
PROCS=$5;
ID=$6;

#location of ngsplot executable
PCODE="/net/isi-scratch/giuseppe/tools/ngsplot/bin";
PRSCRIPT="/net/isi-scratch/giuseppe/tools/R-3.1.0/bin/Rscript";
PFANTOM_ENHANCERS="/net/isi-mirror/FANTOM5/Enhancers/hg19_enhancers.bed";

GENOME="hg19";
COMMON="-P ${PROCS} -H ${SHADING} -D ${PDB} -LWD 1";
#-RB 0.02		


#SCRIPT=ngsplot_fantom_en_${ID}.sh;
#echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
#echo '' >>${PDATA}/${SCRIPT};
#echo "${PRSCRIPT} ${PCODE}/ngs.plot.r -G ${GENOME} -R bed -C ${PDATA}/${CONFIG_FILE} -E ${PFANTOM_ENHANCERS} -O ${PDATA}/${ID}_fantom_en -T ${ID}_fantom_en ${COMMON}" >> ${PDATA}/${SCRIPT};	
#nice -5 qsub -e ${PDATA}/ngsplot_fantom_en_${ID}.err -o ${PDATA}/ngsplot_fantom_en_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
#rm ${PDATA}/${SCRIPT};  

SCRIPT=ngsplot_tss_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "${PRSCRIPT} ${PCODE}/ngs.plot.r -G ${GENOME} -R tss -C ${PDATA}/${CONFIG_FILE} -O ${PDATA}/${ID}_tss -T ${ID}_tss ${COMMON}" >> ${PDATA}/${SCRIPT};
nice -5 qsub -e ${PDATA}/ngsplot_tss_${ID}.err -o ${PDATA}/ngsplot_tss_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};  

SCRIPT=ngsplot_tes_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "${PRSCRIPT} ${PCODE}/ngs.plot.r -G ${GENOME} -R tes -C ${PDATA}/${CONFIG_FILE} -O ${PDATA}/${ID}_tes -T ${ID}_tes ${COMMON}" >> ${PDATA}/${SCRIPT};
nice -5 qsub -e ${PDATA}/ngsplot_tes_${ID}.err -o ${PDATA}/ngsplot_tes_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};  

SCRIPT=ngsplot_enh_pc_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "${PRSCRIPT} ${PCODE}/ngs.plot.r -G ${GENOME} -R enhancer -C ${PDATA}/${CONFIG_FILE} -O ${PDATA}/${ID}_enh_pc -F Gm12878 -T ${ID}_enhancer_pc ${COMMON}" >> ${PDATA}/${SCRIPT};
nice -5 qsub -e ${PDATA}/ngsplot_enh_pc_${ID}.err -o ${PDATA}/ngsplot_enh_pc_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};  

SCRIPT=ngsplot_enh_lnRNA_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "${PRSCRIPT} ${PCODE}/ngs.plot.r -G ${GENOME} -R enhancer -C ${PDATA}/${CONFIG_FILE} -O ${PDATA}/${ID}_enh_lnRNA -F Gm12878,lincRNA -T ${ID}_enhancer_lnRNA ${COMMON}" >> ${PDATA}/${SCRIPT};
nice -5 qsub -e ${PDATA}/ngsplot_enh_lnRNA_${ID}.err -o ${PDATA}/ngsplot_enh_lnRNA_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};  

SCRIPT=ngsplot_dhs_pp_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "${PRSCRIPT} ${PCODE}/ngs.plot.r -G ${GENOME} -R dhs -C ${PDATA}/${CONFIG_FILE} -O ${PDATA}/${ID}_dhs_pp -F Gm12878 -T ${ID}_dhs_pp ${COMMON}" >> ${PDATA}/${SCRIPT};
nice -5 qsub -e ${PDATA}/ngsplot_dhs_pp_${ID}.err -o ${PDATA}/ngsplot_dhs_pp_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};  

SCRIPT=ngsplot_dhs_p3k_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "${PRSCRIPT} ${PCODE}/ngs.plot.r -G ${GENOME} -R dhs -C ${PDATA}/${CONFIG_FILE} -O ${PDATA}/${ID}_dhs_p3k -F Gm12878,Promoter3k -T ${ID}_dhs_p3k  ${COMMON}" >> ${PDATA}/${SCRIPT};
nice -5 qsub -e ${PDATA}/ngsplot_dhs_p3k_${ID}.err -o ${PDATA}/ngsplot_dhs_p3k_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};  

SCRIPT=ngsplot_cgi_pp_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "${PRSCRIPT} ${PCODE}/ngs.plot.r -G ${GENOME} -R cgi -C ${PDATA}/${CONFIG_FILE} -O ${PDATA}/${ID}_cgi_pp -T ${ID}_cgi_pp ${COMMON}" >> ${PDATA}/${SCRIPT};
nice -5 qsub -e ${PDATA}/ngsplot_cgi_pp_${ID}.err -o ${PDATA}/ngsplot_cgi_pp_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};  

SCRIPT=ngsplot_cgi_p3k_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "${PRSCRIPT} ${PCODE}/ngs.plot.r -G ${GENOME} -R cgi -C ${PDATA}/${CONFIG_FILE} -O ${PDATA}/${ID}_cgi_p3k -F Promoter3k -T ${ID}_cgi_p3k ${COMMON}" >> ${PDATA}/${SCRIPT};
nice -5 qsub -e ${PDATA}/ngsplot_cgi_p3k_${ID}.err -o ${PDATA}/ngsplot_cgi_p3k_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};  

SCRIPT=ngsplot_genebody_pc_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "${PRSCRIPT} ${PCODE}/ngs.plot.r -G ${GENOME} -R genebody -C ${PDATA}/${CONFIG_FILE} -O ${PDATA}/${ID}_genebody_pc -T ${ID}_genebody_pc ${COMMON}" >> ${PDATA}/${SCRIPT};
nice -5 qsub -e ${PDATA}/ngsplot_genebody_pc_${ID}.err -o ${PDATA}/ngsplot_genebody_pc_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};  

SCRIPT=ngsplot_genebody_lnRNA_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "${PRSCRIPT} ${PCODE}/ngs.plot.r -G ${GENOME} -R genebody -C ${PDATA}/${CONFIG_FILE} -O ${PDATA}/${ID}_genebody_lnRNA -F lincRNA -T ${ID}_genebody_lnRNA ${COMMON}" >> ${PDATA}/${SCRIPT};
nice -5 qsub -e ${PDATA}/ngsplot_genebody_lnRNA_${ID}.err -o ${PDATA}/ngsplot_genebody_lnRNA_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};  

SCRIPT=ngsplot_genebody_miRNA_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "${PRSCRIPT} ${PCODE}/ngs.plot.r -G ${GENOME} -R genebody -C ${PDATA}/${CONFIG_FILE} -O ${PDATA}/${ID}_genebody_miRNA -F miRNA -T ${ID}_genebody_miRNA ${COMMON}" >> ${PDATA}/${SCRIPT};
nice -5 qsub -e ${PDATA}/ngsplot_genebody_miRNA_${ID}.err -o ${PDATA}/ngsplot_genebody_miRNA_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};
