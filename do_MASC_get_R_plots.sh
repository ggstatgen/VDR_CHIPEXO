#!/bin/bash

#script to automate R plotting of the MASC output.
#It requires a directory of .txt output files from MASC
#It should output a labelled .svg file, one per sample
#Optionally it could output measures similar to what I get from phantompeak tools

#Typical input table
#d       Correlation     Mean Correlation        MaSC    Mean MaSC
#0       0.0143002149653467      0.01445318333667        0.00798208340801875     0.00810406205533825
#1       0.0115085329146738      0.0147361394309213      0.00498761191352638     0.00838730100266387
#2       0.019964151771269       0.0149150190416444      0.0140399983476863      0.00856157431029565
#3       0.0159526231387603      0.0150905583220255      0.00971045920514277     0.00872910413789112
#4       0.0125528832004852      0.015335302887763       0.00602263452801052     0.00897240819618516
#5       0.0124406940294852      0.0152406046670979      0.00588158492964479     0.00885313525038763
#6       0.0164338759964293      0.0154350758533572      0.0100867346866175      0.00902690404529357
#The spaces in the titles will confuse R, remove header?

#I want to plot at least "Correlation" and "MaSC". The means are not too important, maybe I'll plot them only for the POOLED file


if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <PATH> <TRIM> <SVG_W> <SVG_H>"
        echo "<PATH> data path"
	echo "<TRIM> d max you want to visualise in the plot (eg 100)"
	echo "<SVG_W> width of the svg output file  - eg 15"
	echo "<SVG_H> height of the svg output file - eg  8"
        exit
fi

PDATA=$1;
TRIM=$2;
SVG_W=$3;
SVG_H=$4;
PCODE="/net/isi-scratch/giuseppe/tools/R-3.1.0/bin/";

for FILE in ${PDATA}/*.txt;
	do ID=`basename ${FILE} ".txt"`;
	SCRIPT=${ID}.r;
	#BASENAME=MaSC_${ID};

	echo "library(ggplot2)" >> ${PDATA}/${SCRIPT};
	echo "data <- read.csv('${FILE}',sep= \"\t\", header =  T,check.names = F)" >> ${PDATA}/${SCRIPT};
	#stats about d_cl and correlation vals
	#echo 'min_cc <- abs(min(abs(data$Correlation)))' >> ${PDATA}/${SCRIPT};
	#echo 'min_cc_masc <- abs(min(abs(data$MaSC)))' >> ${PDATA}/${SCRIPT};
        echo 'min_cc <- min(data$Correlation)' >> ${PDATA}/${SCRIPT};
        echo 'min_cc_masc <- min(data$MaSC)' >> ${PDATA}/${SCRIPT};

	#for the chipseq: what is the estmated d for the phantom peak?
	#the peak is small, I need to trim soonafter the readlength (36 or 40)
	echo 'data_pp <- data[35:45,]' >> ${PDATA}/${SCRIPT};
        echo 'd_pp <- data_pp$d[which(data_pp$Correlation == max(data_pp$Correlation))]' >> ${PDATA}/${SCRIPT};

	#chip exo: trim the data up to the phantompeak for the following?
	#echo 'data_tr <- data[0:35,]' >> ${PDATA}/${SCRIPT};

	#mappability-corrected NSC
	echo 'd_cs_masc <- data$d[which(data$MaSC == max(data$MaSC))]' >> ${PDATA}/${SCRIPT};
	echo 'cc_d_cs_masc <- max(data$MaSC)' >> ${PDATA}/${SCRIPT};
	echo 'NSC_masc <- cc_d_cs_masc / min_cc_masc' >> ${PDATA}/${SCRIPT};

	#naive NSC
        echo 'd_cs <- data$d[which(data$Correlation == max(data$Correlation))]' >> ${PDATA}/${SCRIPT};
	echo 'cc_d_cs <- max(data$Correlation)' >> ${PDATA}/${SCRIPT};
        echo 'NSC <- cc_d_cs / min_cc' >> ${PDATA}/${SCRIPT};

	#avgs
        echo 'd_cs_avg <- data$d[which(data$`Mean Correlation` == max(data$`Mean Correlation`))]' >> ${PDATA}/${SCRIPT};
        echo 'd_cs_masc_avg <- data$d[which(data$`Mean MaSC` == max(data$`Mean MaSC`))]' >> ${PDATA}/${SCRIPT};
	
	#print this stuff
	echo "fileConn<-file('${PDATA}/${ID}.Rdata')" >> ${PDATA}/${SCRIPT};
	#echo 'writeLines(c(paste(d_cs,d_cs_masc,min_cc,min_cc_masc,NSC,NSC_masc)), fileConn)' >> ${PDATA}/${SCRIPT};
	#echo 'writeLines(c(paste(d_cs,d_cs_masc,min_cc_masc,NSC_masc,)), fileConn)' >> ${PDATA}/${SCRIPT};
	#echo  'writeLines(c(paste(d_cs_masc,cc_d_cs_masc, min_cc_masc,NSC_masc)), fileConn)' >> ${PDATA}/${SCRIPT};
	echo 'writeLines(c(paste(d_pp,d_cs,d_cs_masc,NSC,NSC_masc)), fileConn)' >> ${PDATA}/${SCRIPT};
	echo "close(fileConn)" >> ${PDATA}/${SCRIPT};

	#trim data to the first TRIM for plotting
	echo "data_tr <- data[0:$TRIM,]" >> ${PDATA}/${SCRIPT};

	#stack data in one column
	#MEAN COLUMNS = NO:
	#echo 'data2 <- cbind(data[gl(nrow(data),1,2*nrow(data)),1],stack(data[,c(2,4)]))' >> ${PDATA}/${SCRIPT};
	#MEAN COLUMNS = YES:
	#do an if here to check if there are less columns than TRIM
	#echo "data2 <- cbind(data[gl(nrow(data),1,4*nrow(data)),1],stack(data[,2:5]))" >> ${PDATA}/${SCRIPT};
	echo "data2 <- cbind(data_tr[gl(nrow(data_tr),1,4*nrow(data_tr)),1],stack(data_tr[,2:5]))" >> ${PDATA}/${SCRIPT};
	
	#print
	echo 'names(data2)[1] <- "numero"' >> ${PDATA}/${SCRIPT};
	
	#MEAN COLUMNS = NO:
	#echo 'ggplot(data=data2, aes(x=numero, y=values, group = ind, colour = ind)) + 
        #geom_line() + 
        #labs(x = "Strand Shift (bp)", y = "Cross-correlation") +
        #scale_color_discrete(name = "Method", labels = c("Mpb Corr", "Naive")) +
	#geom_vline(xintercept=c(d_cs_masc,6,9,12), linetype="dotted") +
        #theme_bw(base_size = 12, base_family = "Helvetica")' >> ${PDATA}/${SCRIPT};

	#MEAN COLUMNS = YES:
	#        ylim(0,0.002) +
	#scale_y_log10(limits = c(0.0001,0.006)) +
	echo 'ggplot(data=data2, aes(x=numero, y=values, group = ind, colour = ind)) + 
        geom_line() + 
        labs(x = "Strand Shift (bp)", y = "Cross-correlation") +
        scale_color_discrete(name = "Method", labels = c("Naive", "Mpb Corr",  "Avg Naive", "Avg Mpb Corr")) +
        geom_vline(xintercept=c(40), linetype="dotted") +
        theme_bw(base_size = 15, base_family = "Helvetica")' >> ${PDATA}/${SCRIPT};
	
	echo "ggsave(file = '${PDATA}/${ID}.svg', width=${SVG_W}, height=${SVG_H})" >> ${PDATA}/${SCRIPT}; 

	#execute script
	#${PCODE}/Rscript ${PDATA}/${SCRIPT};
	#rm ${PDATA}/${SCRIPT};
done
