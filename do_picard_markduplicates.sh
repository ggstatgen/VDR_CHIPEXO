#!/bin/bash

#Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. All records are then written to the output file with the duplicate records flagged.

#ATTENTION: if you sorted your bams using samtools, and try this program on them, you'll get an exception: "The file xxx is not coordinate sorted". This is because Picard expects a flag
#to be set to indicate sorted status. Samtools does not set it. To make this go on, use the AS=true flag


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input .bam files"
	exit
fi

PDATA=$1;

#------------------------
#hard coded paths
#------------------------
PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.117";
TMP="/tmp";

#              REMOVE_DUPLICATES=TRUE \\
for FILE in ${PDATA}/*.bam;
        do
        FILE_ID=`basename ${FILE} ".bam"`;

        SCRIPT=picard_md_${FILE_ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        echo "java -Xmx4g -Djava.io.tmpdir=${TMP} \\
              -jar ${PCODE_PICARD}/MarkDuplicates.jar \\
              INPUT=${FILE} \\
              OUTPUT=${PDATA}/${FILE_ID}_md.bam \\
              REMOVE_DUPLICATES=TRUE \\
              AS=true \\
              METRICS_FILE=${PDATA}/picard_md_${FILE_ID}.metrics" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/picard_md_${FILE_ID}.stderr -o ${PDATA}/picard_md_${FILE_ID}.stdout -v "BASH_ENV=~/.bashrc" -q fgu217.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
