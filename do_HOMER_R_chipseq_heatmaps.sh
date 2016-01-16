#!/bin/bash

#13/10/2013
#use HOMER and R to obtain heatmaps of the peaks
#http://biowhat.ucsd.edu/homer/ngs/quantification.html
#INPUTS
#1. peak summit file or peak bed file produced by MACS2 or GEM (tested with MACS2)
#2. bam alignment obtained via bwa or bowtie
#OUTPUTS
#pdf plots with heatmaps like these
#http://crazyhottommy.blogspot.co.uk/2013/04/how-to-make-tss-plot-using-rna-seq-and.html

#INPUT PARAMETERS
#number of best peaks to consider
#homer parameters (size, sampling, etc)

#steps
#1 select the n best peaks using 
#cat MACS2_picard_merged.sorted_peaks.narrowPeak | sort -k8nr | head -n 1000 
#2 create homer bed peaks using 
#do_npk_to_homer.pl
#3 HOMER create tag directory with homer
#4 HOMER annotate peaks
#5 vim outfile and leave only data
#6 run R script on the clean HOMER output 

#makeTagDirectory
#Usage: makeTagDirectory <directory> <alignment file 1> [file 2] ... [options]

#I assume the file naming convention is
#BED: MACS2_[ID]_summits.bed
#or
#BED: MACS2_[ID]_peaks.narrowPeak
#and
#BAM: [ID].bam

#example 1
#BED1: MACS2_VDR_GM06986_reads_summits.bed
#BED2: MACS2_VDR_GM06986_reads_peaks.narrowPeak
#BAM: VDR_GM06986_reads.bam
#example 2
#BED:
#BAM:

if [ ! $# == 6 ]; then
        echo "Usage: `basename $0` <PATH_INTERVAL> <INTERVAL_TYPE> <PATH_BAM> <HOMER_SIZE> <HOMER_HIST_BIN> <HOMER_LENGTH>"
        echo "<PATH_INTERVAL> path to the directory containing the interval files (either MACS2 summits or MACS2 narrowPeak)"
	echo "<INTERVAL_TYPE> [n|s] (narrowpeak (n) or summit(s))"
	echo "<PATH_BAM> path to the directory containing the bam aligments"
	echo "<HOMER_SIZE> (e.g. 1000)"
	echo "<HOMER_HIST_BIN> (e.g. 2)"
	echo "<HOMER_LENGTH> (e.g. 0)"
        exit
fi

P_BED=$1;
INT_TYPE=$2;
P_BAM=$3;
HOMER_SIZE=$4;
HOMER_HIST=$5;
HOMER_LEN=$6;

PREFIX="MACS2";
POSTFIX_NPK="peaks.narrowPeak";
POSTFIX_SUMMIT="summits.bed";
TAGDIR=${P_BAM}/"d_HOMER_TAGS";

PATH_R="/net/isi-scratch/giuseppe/tools/R-3.0.1/bin";
PATH_HOMER="/net/isi-scratch/giuseppe/tools/HOMER/bin";
PATH_SCRIPTS="/net/isi-backup/giuseppe/scripts";
GENOME_SIZE="2540757438";

for BAM_FILE in ${P_BAM}/*.bam;
        do 
	ID=`basename ${BAM_FILE} ".bam"`;
        SCRIPT=HOMER_ap_${ID}.sh;
	
	#build the interval file name based on option 2 and on the given bam ID
	if [[ ${INT_TYPE} -eq n ]];
	then
        	INTERVAL_FILE=${P_BED}/${PREFIX}_${ID}_${POSTFIX_NPK};
	elif [[ ${INT_TYPE} -eq s ]];
	then
        	INTERVAL_FILE=${P_BED}/${PREFIX}_${ID}_${POSTFIX_SUMMIT};
	else
		echo '<INTERVAL_TYPE> field not recognised. Aborting.';
		exit;
	fi

        echo '#!/bin/bash' >>${P_BED}/${SCRIPT};
        echo '' >>${P_BED}/${SCRIPT};
	#------------------------------------
	#1. create HOMER-compatible bed files
	#------------------------------------
	echo "${PATH_SCRIPTS}/do_npk_to_HOMER_bed.pl -i=${INTERVAL_FILE} > ${P_BED}/HOMER_${ID}.intervals.bed" >>${P_BED}/${SCRIPT}; #should also work with summit file
	#----------------------------------------------------
	#2. create HOMER-compratible tag directory from bam file - only if the directory does not exist
	#----------------------------------------------------
	if [ ! -d "$TAGDIR" ]; then
		# Control will enter here if $DIRECTORY doesn't exist.
		echo "${PATH_HOMER}/makeTagDirectory ${TAGDIR} ${BAM_FILE} -genome hg19 > ${TAGDIR}/HOMER_${ID}_makeTagDirectory.log" >>${P_BED}/${SCRIPT};
	fi
 	#------------------------
	#3. HOMER - annotatepeaks
	#-----------------------
	#removed:               -len ${HOMER_LEN} \\
	echo "${PATH_HOMER}/annotatePeaks.pl ${P_BED}/HOMER_${ID}.intervals.bed  hg19 \\
              -gsize ${GENOME_SIZE} \\
              -d ${TAGDIR} \\
              -hist ${HOMER_HIST} \\
              -ghist \\
              -noadj \\
              -len ${HOMER_LEN} \\
              -size ${HOMER_SIZE} > ${P_BED}/HOMER_${ID}.matrix" >>${P_BED}/${SCRIPT};
	nice -5 qsub -e ${P_BED}/HOMER_${ID}.err -q newnodes.q ${P_BED}/${SCRIPT};
        rm ${P_BED}/${SCRIPT};
done

#check that the former set of jobs has finished
#sleep if there are still jobs executing
JOB_STATUS=$(qstat);
if [[ $JOB_STATUS =~ HOMER_ap ]] ;
then
        COMPLETE=0;
else
        COMPLETE=1;
fi

while [ "${COMPLETE}" = "0" ]
do
        echo "not all HOMER_ap jobs complete. Waiting 10 seconds.."
        sleep 10;

        JOB_STATUS=$(qstat);
        if [[ $JOB_STATUS =~ HOMER_ap ]] ;
        then
                COMPLETE=0;
        else
                COMPLETE=1;
        fi
done
echo "All jobs completed, continuing to R analysis..";
echo " ";


#for the full script look at R_chipseq_heatmaps.R
for FILE in ${P_BED}/*.matrix;
	do ID=`basename ${FILE} ".matrix"`;
	SCRIPT=HOMER_r_${ID}.r;

	echo 'library(ggplot2)' >> ${P_BED}/${SCRIPT};
	echo 'library(pheatmap)' >> ${P_BED}/${SCRIPT};
	echo "d1 <- read.table('${FILE}', header=T)" >> ${P_BED}/${SCRIPT};
	# heatmap.2 works only on matrix, turn the dataframe to matrix, and add the TSS id as the row name
	echo 'm1<- as.matrix( d1[,2:ncol(d1)])' >> ${P_BED}/${SCRIPT};
	echo 'rownames(m1)<- d1$Gene' >> ${P_BED}/${SCRIPT};
	# log2 transform the raw counts
	echo 'm1<- log2(m1+1)' >> ${P_BED}/${SCRIPT};
	#reorder the data from largest to narrowest peak
	#save sum of counts in last column
	echo 'm.row.sum <- cbind(m1, rowSums(m1))' >> ${P_BED}/${SCRIPT};
	echo 'o1<- rev(order(m.row.sum[,ncol(m.row.sum)]))' >> ${P_BED}/${SCRIPT};
	echo 'm.row.sum<- m.row.sum[o1,]' >> ${P_BED}/${SCRIPT};
	#colour palette white to red
	echo 'bk = unique(c(seq(-0.1,3, length=100),seq(3,10.35,length=100)))' >> ${P_BED}/${SCRIPT};
	echo 'hmcols<- colorRampPalette(c("white","red"))(length(bk)-1)' >> ${P_BED}/${SCRIPT};

	echo "png('${P_BED}/R_${ID}.heatmap.png', width=600, height = 1600)" >> ${P_BED}/${SCRIPT};
	#echo "pdf('${P_BED}/HOMER_R_${ID}_heatmap.pdf')" >> ${P_BED}/${SCRIPT};
	echo 'pheatmap( m.row.sum[,1:(ncol(m.row.sum)-1)], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)' >> ${P_BED}/${SCRIPT};
	echo 'dev.off()' >> ${P_BED}/${SCRIPT};

	#execute script
	${PATH_R}/Rscript ${P_BED}/${SCRIPT};
done
