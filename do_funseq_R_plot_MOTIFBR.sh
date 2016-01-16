#!/bin/bash

#Script to automate R ggplot2 plotting of MOTIFBR funseq2 results
#Inputs are a list of Rdata files produced by do_funseq_collect_MOTIFBR_gtph_cluster.sh

#typical steps
#> library(ggplot2)
#> library(reshape2)
#> data_CREB1_MA0018.1 <- read.table("Output_noDBRECUR_CREB1_MA0018.1_concordance_all.Rdata",sep="\t",header=T)
#> data_CREB1_MA0018.1_m <- melt(data_CREB1_MA0018.1)
#data_CREB1_MA0018.1_m$MOTIF_POS <- factor(data_CREB1_MA0018.1_m$MOTIF_POS, levels=data_CREB1_MA0018.1$MOTIF_POS)
#> ggplot(data_CREB1_MA0018.1_m, aes(x=MOTIF_POS,y=value,fill=factor(variable))) + geom_bar(stat="identity",width=.5) +scale_y_continuous(limits=c(0, 20),breaks=seq(from=0,to=20,by=2)) + theme_minimal() + scale_fill_manual(values=c("#FF3300", "#999999"))
#> ggsave("data_CREB1_MA0018.1.svg", width=15)
#> ggsave("data_CREB1_MA0018.1.png", width=15)


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "<PATH> data path for the .Rdata files";
	echo "NOTE: input files are currently expected to be in the form Output_noDBRECUR_<JASPARTFNAME>_<JASPARTFID>_concordance_all.Rdata";
	echo "with the two jaspar field at the 3,4 position. Modify script if not";
        exit
fi

PDATA=$1;
PCODE="/net/isi-scratch/giuseppe/tools/R-3.1.0/bin";

#Output_noDBRECUR_Arid3a_MA0151.1_concordance_all.Rdata
for FILE in ${PDATA}/*.Rdata;
        do
        SEED=`echo ${FILE} | grep -Po "RECUR_(.*)_concordance"`;
        TFBASE=`echo ${SEED} | cut -d _ -f 2`;
	TF=`echo ${TFBASE/-/}`; #remove any hyphens from the TF name, because R does not like variable names with hyphens
        
	ID=`echo ${SEED} | cut -d _ -f 3`;
	SCRIPT=plot_MOTIFBR_${TF}_${ID}.r
	
	echo " ";
	echo "Processing $TF - $ID:";
	echo " ";

	echo "library(ggplot2)"  >> ${PDATA}/${SCRIPT};
	echo "library(reshape2)" >> ${PDATA}/${SCRIPT};
	echo "data_${TF}_${ID} <- read.table(\"${FILE}\",sep=\"\t\",header=T)" >> ${PDATA}/${SCRIPT};
	echo "data_${TF}_${ID}_m <- melt(data_${TF}_${ID})" >> ${PDATA}/${SCRIPT};
	echo "data_${TF}_${ID}_m\$MOTIF_POS <- factor(data_${TF}_${ID}_m\$MOTIF_POS, levels=data_${TF}_${ID}\$MOTIF_POS)" >> ${PDATA}/${SCRIPT};
	echo "ggplot(data_${TF}_${ID}_m, aes(x=MOTIF_POS,y=value,fill=factor(variable))) + geom_bar(stat=\"identity\",width=.7) +scale_y_continuous(limits=c(0, 20),breaks=seq(from=0,to=20,by=2)) + theme_minimal() + scale_fill_manual(values=c(\"#FF3300\", \"#999999\"))" >> ${PDATA}/${SCRIPT};
	echo "ggsave(\"${PDATA}/data_${TF}_${ID}.svg\", width=15)" >> ${PDATA}/${SCRIPT};
	echo "ggsave(\"${PDATA}/data_${TF}_${ID}.png\", width=15)" >> ${PDATA}/${SCRIPT};

	#execute script
	${PCODE}/Rscript ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT};
done
