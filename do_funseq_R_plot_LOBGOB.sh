#!/bin/bash

#Script to automate R ggplot2 plotting of LOBGOB  funseq2 results
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


if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <YSCALE>"
        echo "<PATH> data path for the .Rdata files";
	echo "<YSCALE> y scale to try - you will need to set this dynamically in R..for now try 100";
	echo "NOTE: input files are currently expected to be in the form Output_noDBRECUR_<JASPARTFNAME>_<JASPARTFID>_phdisruption.Rdata";
	echo "with the two jaspar field at the 3,4 position. Modify script if not";
        exit
fi

PDATA=$1;
YSCALE=$2;
PCODE="/net/isi-scratch/giuseppe/tools/R-3.1.0/bin";

#Output_noDBRECUR_Arid3a_MA0151.1_concordance_all.Rdata
for FILE in ${PDATA}/*.Rdata;
        do
        SEED=`echo ${FILE} | grep -Po "RECUR_(.*)_phdisruption"`;
        TFBASE=`echo ${SEED} | cut -d _ -f 2`;
	TF=`echo ${TFBASE/-/}`; #remove any hyphens from the TF name, because R does not like variable names with hyphens
        
	ID=`echo ${SEED} | cut -d _ -f 3`;
	SCRIPT=plot_MOTIFBR_${TF}_${ID}.r
	
	echo " ";
	echo "Processing $TF - $ID:";
	echo " ";

	echo "library(ggplot2)"  >> ${PDATA}/${SCRIPT};
	echo "data_${TF}_${ID} <- read.table(\"${FILE}\",sep=\"\t\",header=T)" >> ${PDATA}/${SCRIPT};
	#subset(data,!duplicated(data$ID))
	echo "data_${TF}_${ID}_dedup <- subset(data_${TF}_${ID},!duplicated(data_${TF}_${ID}\$position))" >> ${PDATA}/${SCRIPT};
	echo "data_${TF}_${ID}\$position <- factor(data_${TF}_${ID}\$position, levels=data_${TF}_${ID}_dedup\$position)" >> ${PDATA}/${SCRIPT};
	
	#if you want again the gray for the SYM, it's
	echo "ggplot(data_${TF}_${ID}, aes(y=((ref+1)/(alt+1)), x=position))  + geom_jitter(size=5, alpha = I(.7), position = position_jitter(width = .2,height=0), aes(colour = type)) + theme_minimal() + ggtitle(\"Impact of LOB VDR-BVs on binding affinity at ${TF}(${ID}) consensus motif\") + ylab(\"rc(anc)+1/rc(alt)+1\") +  xlab(\"Position in Consensus Motif\") +  scale_y_continuous(limits = c(1, ${YSCALE})) +  scale_colour_manual(values=c(\"#3399FF\", \"orange\", \"#FF3300\", \"#0000ff\",\"#999999\"))"  >> ${PDATA}/${SCRIPT};
	echo "ggsave(\"${PDATA}/data_LOB_${TF}_${ID}.svg\", width=15)" >> ${PDATA}/${SCRIPT};
	echo "ggsave(\"${PDATA}/data_LOB_${TF}_${ID}.png\", width=15)" >> ${PDATA}/${SCRIPT};
	
	echo "ggplot(data_${TF}_${ID}, aes(y=((alt+1)/(ref+1)), x=position))  + geom_jitter(size=5, alpha = I(.7),position =  position_jitter(width = .2,height=0), aes(colour = type)) + theme_minimal() + ggtitle (\"Impact of GOB VDR-BVs on binding affinity at ${TF}(${ID}) consensus motif\") + ylab(\"rc(alt)+1/rc(anc)+1\") +  xlab(\"Position in Consensus Motif\") +  scale_y_continuous(limits = c(1, ${YSCALE})) +  scale_colour_manual(values=c(\"#3399FF\", \"orange\", \"#FF3300\", \"#0000ff\",\"#999999\"))"  >> ${PDATA}/${SCRIPT}; 
	echo "ggsave(\"${PDATA}/data_GOB_${TF}_${ID}.svg\", width=15)" >> ${PDATA}/${SCRIPT};
	echo "ggsave(\"${PDATA}/data_GOB_${TF}_${ID}.png\", width=15)" >> ${PDATA}/${SCRIPT};
	
	#execute script
	${PCODE}/Rscript ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT};
done
