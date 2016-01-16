#!/bin/bash

#Takes an input BAM and reference sequence and runs one or more Picard metrics modules at the same time to cut down on I/O. Currently all programs are run with default options and fixed output extesions, but this may become more flexible in future.

#ATTENTION: if you sorted your bams using samtools, and try this program on them, you'll get an exception: "The file xxx is not coordinate sorted". This is because Picard expects a flag
#to be set to indicate sorted status. Samtools does not set it. To make this go on, use the AS=true flag


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input .bam files"
	echo "Assuming hg19 reference"
	exit
fi

PDATA=$1;

#------------------------
#hard coded paths
#------------------------
PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.117";
PREFERENCE="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";
TMP="/tmp";

for FILE in ${PDATA}/*.bam;
        do
        FILE_ID=`basename ${FILE} ".bam"`;

        SCRIPT=picard_metrics_${FILE_ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        echo "java -Xmx4g -Djava.io.tmpdir=${TMP} \\
              -jar ${PCODE_PICARD}/CollectMultipleMetrics.jar.jar \\
              INPUT=${FILE} \\
              REFERENCE_SEQUENCE=${PREFERENCE} \\
              OUTPUT=${PDATA}/${FILE_ID}_picard.metrics \\
              PROGRAM=CollectAlignmentSummaryMetrics \\
              PROGRAM=QualityScoreDistribution" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/picard_metrics_${FILE_ID}.stderr -o ${PDATA}/picard_metrics_${FILE_ID}.stdout -v "BASH_ENV=~/.bashrc" -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
