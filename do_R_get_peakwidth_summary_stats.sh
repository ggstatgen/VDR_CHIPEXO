#!/bin/bash

#This gives me a file with per-sample summary statistics of peak width
#It expects a directory with .hist files generated with do_bed_get_peak_sizes_cluster.sh 

#input is 
#GM07029 24
#GM07029 109
#GM07029 59
#GM07029 49
#....

#it basically calls r summary() and returns, for each sample, a line: 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  12.00   12.00   12.00   24.71   28.00  345.00

#you can gather these and make a table

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "<PATH> data path for the .hist files"
        exit
fi

PDATA=$1;
PCODE="/net/isi-scratch/giuseppe/tools/R-3.1.0/bin/";

for FILE in ${PDATA}/*.hist;
	do ID=`basename ${FILE} ".hist"`;
	SCRIPT=summarystats_${ID}.r;
	#BASENAME=MaSC_${ID};

	echo "data <- read.table('${FILE}')" >> ${PDATA}/${SCRIPT};
	#echo 'p1 <- summary(data$V1)' >> ${PDATA}/${SCRIPT};
	echo 'p <- summary(data$V2)' >>  ${PDATA}/${SCRIPT};
	#print this stuff

	echo "write.csv(t(as.matrix(p)), file='${PDATA}/${ID}.Rdata', quote=FALSE, row.names=F)" >> ${PDATA}/${SCRIPT};

        #echo "fileConn<-file('${PDATA}/${ID}.Rdata')" >> ${PDATA}/${SCRIPT};
        #echo 'writeLines(c(paste(p)), fileConn)' >> ${PDATA}/${SCRIPT};
        #echo "close(fileConn)" >> ${PDATA}/${SCRIPT};

        #execute script
        ${PCODE}/Rscript ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done

echo "ID:Min.,1st Qu.,Median,Mean,3rd Qu.,Max." > ${PDATA}/peak_width_table.csv
grep -v "Min" ${PDATA}/*.Rdata >> ${PDATA}/peak_width_table.csv

