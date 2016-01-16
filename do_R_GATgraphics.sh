#!/bin/bash

#script to automate R plotting of the GAT output.
#It requires a .tsv file output by GAT
#it outputs a labelled pdf with a horizontal bar plot and pie chart


if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <Y_AXIS> <PLOT_TITLE>"
        echo "<PATH> data path"
	echo "<Y_AXIS> one of [l2fold|fold]"
	echo "<PLOT_TITLE> generic plot title"
        exit
fi

PDATA=$1;
YAXIS=$2;
PLOT_TITLE=$3;
PCODE="/net/isi-scratch/giuseppe/tools/R-3.1.0/bin";

for FILE in ${PDATA}/*.tsv;
        do ID=`basename ${FILE} ".tsv"`;
        SCRIPT=${ID}.r;

        echo "library(ggplot2)" >> ${PDATA}/${SCRIPT};
        echo "myData<-read.table('${FILE}', header=T)" >> ${PDATA}/${SCRIPT};
        echo 'ggplot(myData, aes(y='${YAXIS}', x=reorder(annotation,'${YAXIS}')), fill=annotation) +
        geom_bar(stat="identity") +
        coord_flip() +
        ggtitle('\"${PLOT_TITLE}\"') +
        xlab("Annotation (Ensembl 72)") +
	theme_bw(base_size = 12, base_family = "Helvetica") -> p' >> ${PDATA}/${SCRIPT};
        echo "pdf('${PDATA}/${ID}_enrichment.pdf')" >> ${PDATA}/${SCRIPT};
        echo "print(p)" >> ${PDATA}/${SCRIPT};
        echo "dev.off()" >> ${PDATA}/${SCRIPT};

        #pie chart
        echo 'slices <- myData$observed' >> ${PDATA}/${SCRIPT};
        echo 'lbls <- myData$annotation' >> ${PDATA}/${SCRIPT};
        echo 'pct <- round(slices/sum(slices)*100)' >> ${PDATA}/${SCRIPT};
        echo 'lbls <- paste(lbls, pct)' >> ${PDATA}/${SCRIPT};
        echo 'lbls <- paste(lbls,"%",sep="")' >> ${PDATA}/${SCRIPT};
        echo "pdf('${PDATA}/${ID}_piechart.pdf')" >> ${PDATA}/${SCRIPT};
        echo 'pie(slices, labels = lbls, main="GAT VDR ChIP-exo")' >> ${PDATA}/${SCRIPT};
        echo "dev.off()" >> ${PDATA}/${SCRIPT};

        #execute script
        ${PCODE}/Rscript ${PDATA}/${SCRIPT};
done



#for FILE in ${PDATA}/*.tsv;
#	do ID=`basename ${FILE} ".tsv"`;
#	SCRIPT=${ID}.r;
#
#	echo "library(ggplot2)" >> ${PDATA}/${SCRIPT};
#	echo "myData<-read.table('${FILE}', header=T)" >> ${PDATA}/${SCRIPT};
#	echo 'ggplot(myData, aes(y=l2fold, x=reorder(annotation,l2fold)), fill=annotation) +
#	geom_bar(stat="identity",colour="black",fill="lightgreen") +
#	coord_flip() +
#	ggtitle('\"${PLOT_TITLE}\"') +
#	ylab("log2(fold change)") +
#	xlab("Annotation (Ensembl 72)") +
#	theme(plot.title = element_text(face = "bold", size=15)) +
#	theme(axis.text.y = element_text(family = "sans", face = "bold", size = 12)) +
#	theme(axis.text.x = element_text(family = "sans", face = "bold", size = 12)) -> p' >> ${PDATA}/${SCRIPT};
#	echo "pdf('${PDATA}/${ID}_enrichment.pdf')" >> ${PDATA}/${SCRIPT};
#	echo "print(p)" >> ${PDATA}/${SCRIPT};
#	echo "dev.off()" >> ${PDATA}/${SCRIPT};
#
#	#pie chart
#	echo 'slices <- myData$observed' >> ${PDATA}/${SCRIPT};
#	echo 'lbls <- myData$annotation' >> ${PDATA}/${SCRIPT};
#	echo 'pct <- round(slices/sum(slices)*100)' >> ${PDATA}/${SCRIPT};
#	echo 'lbls <- paste(lbls, pct)' >> ${PDATA}/${SCRIPT};
#	echo 'lbls <- paste(lbls,"%",sep="")' >> ${PDATA}/${SCRIPT}; 
#	echo "pdf('${PDATA}/${ID}_piechart.pdf')" >> ${PDATA}/${SCRIPT};
#	echo 'pie(slices, labels = lbls, main="GAT VDR ChIP-exo")' >> ${PDATA}/${SCRIPT};
#	echo "dev.off()" >> ${PDATA}/${SCRIPT};
#
	#execute script
#	${PCODE}/Rscript ${PDATA}/${SCRIPT};
#done
