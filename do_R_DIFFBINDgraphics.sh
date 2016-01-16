#!/bin/bash

#script to automate R Diffbind analysis
#It requires a .csv table with the DiffBind data and a few command line arguments
#it should output a number of plots in pdf

#typical steps
#library(DiffBind)
#vdr = dba(sampleSheet=csvfile)
#olap.rate = dba.overlap(vdr, mode=DBA_OLAP_RATE)
#plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')
#maybe a few of these
#vdr_o2 = dba.count(vdr, minOverlap=2, bParallel=TRUE, score=DBA_SCORE_TMM_READS_EFFECTIVE,bRemoveDuplicates=FALSE)
#vdr_o3 = dba.count(vdr, minOverlap=2, bParallel=TRUE, score=DBA_SCORE_TMM_READS_EFFECTIVE,bRemoveDuplicates=FALSE)
#vdr_o2_nodup = dba.count(vdr, minOverlap=2, bParallel=TRUE, score=DBA_SCORE_TMM_READS_EFFECTIVE,bRemoveDuplicates=TRUE)
#vdr = dba.contrast(vdr, categories=DBA_CONDITION)
#vdr = dba.contrast(vdr, categories=DBA_CONDITION)
#vdr = dba.contrast(vdr, categories=DBA_CONDITION)
#vdr = dba.contrast(vdr, categories=DBA_CONDITION)
#vdr_o1_full = dba.analyze(vdr_o1_full,method=c(DBA_EDGER,DBA_DESEQ), bFullLibrarySize=TRUE)
#dba.plotPCA(vdr,contrast=1,th=.05)
#dba.plotPCA(vdr_o1,contrast=2,method=DBA_EDGER,th=.0001, attributes=c(DBA_FACTOR,DBA_CONDITION), b3D=TRUE)

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <OVERLAP>"
        echo "<PATH> data path"
	echo "<OVERLAP> minimum overlap for peaks across samples"
        exit
fi

PDATA=$1;
MIN_OLP=$2;

PCODE="/net/isi-scratch/giuseppe/tools/R-3.0.1/bin";

for FILE in ${PDATA}/*.csv;
	do ID=`basename ${FILE} ".csv"`;
	SCRIPT=${ID}.r;

	echo "library(DiffBind)" >> ${PDATA}/${SCRIPT};
	echo "vdr = dba(sampleSheet='${FILE}')" >> ${PDATA}/${SCRIPT};
	echo "olap.rate = dba.overlap(vdr, mode=DBA_OLAP_RATE)" >> ${PDATA}/${SCRIPT};
	#output olap plot to pdf
	echo "pdf('${PDATA}/diffbind_${ID}_olap.pdf')" >> ${PDATA}/${SCRIPT};
	echo "plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')" >> ${PDATA}/${SCRIPT};
	echo "dev.off()" >> ${PDATA}/${SCRIPT};
	#counts
	echo 'vdr_count = dba.count(vdr, minOverlap='${MIN_OLP}', bParallel=TRUE, score=DBA_SCORE_TMM_READS_EFFECTIVE,bRemoveDuplicates=FALSE)'  >> ${PDATA}/${SCRIPT};
	echo 'vdr_count_nodup = dba.count(vdr, minOverlap='${MIN_OLP}', bParallel=TRUE, score=DBA_SCORE_TMM_READS_EFFECTIVE,bRemoveDuplicates=TRUE)' >> ${PDATA}/${SCRIPT};
	#set contrast
	echo "vdr_count = dba.contrast(vdr_count, categories=DBA_CONDITION)" >> ${PDATA}/${SCRIPT}; #maybe allow to change category?
	echo "vdr_count_nodup = dba.contrast(vdr_count_nodup, categories=DBA_CONDITION)" >> ${PDATA}/${SCRIPT};
	#diff analysis
	#add #bTagwise=FALSE if there are less than 3>samples per contrast
	echo 'vdr_count = dba.analyze(vdr_count,method=c(DBA_EDGER,DBA_DESEQ))' >> ${PDATA}/${SCRIPT};
	echo 'vdr_count_nodup = dba.analyze(vdr_count_nodup,method=c(DBA_EDGER,DBA_DESEQ))' >> ${PDATA}/${SCRIPT};
	#get correlation plots (8 pdf files)
	echo "pdf('${PDATA}/diffbind_${ID}_analysis_heatmaps.pdf')" >>  ${PDATA}/${SCRIPT};
	echo "dba.plotHeatmap(vdr_count, contrast = 1, th=.05, method=DBA_EDGER)" >> ${PDATA}/${SCRIPT};
	echo "dba.plotHeatmap(vdr_count, contrast = 1, th=.01, method=DBA_EDGER)" >> ${PDATA}/${SCRIPT};
	echo "dba.plotHeatmap(vdr_count, contrast = 1, th=.05, method=DBA_DESEQ)" >> ${PDATA}/${SCRIPT};
	echo "dba.plotHeatmap(vdr_count, contrast = 1, th=.01, method=DBA_DESEQ)" >> ${PDATA}/${SCRIPT};
	#same for the nodup analysis
	echo "dba.plotHeatmap(vdr_count_nodup, contrast = 1, th=.05, method=DBA_EDGER)" >> ${PDATA}/${SCRIPT};
	echo "dba.plotHeatmap(vdr_count_nodup, contrast = 1, th=.01, method=DBA_EDGER)" >> ${PDATA}/${SCRIPT};
	echo "dba.plotHeatmap(vdr_count_nodup, contrast = 1, th=.05, method=DBA_DESEQ)" >> ${PDATA}/${SCRIPT};
	echo "dba.plotHeatmap(vdr_count_nodup, contrast = 1, th=.01, method=DBA_DESEQ)" >> ${PDATA}/${SCRIPT};
	echo "dev.off()" >> ${PDATA}/${SCRIPT};

	#execute script
	${PCODE}/Rscript ${PDATA}/${SCRIPT};
done
