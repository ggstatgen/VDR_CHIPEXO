#!/bin/bash
#script to perform bam>fastq conversion
#bam>sam (samtools) and sam>fastq (picard)
if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input .bam files"
        exit
fi

#VALIDATION_STRINGENCY=LENIENT

PDATA=$1;
PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.102";
#temp
TMP="/tmp";

for FILE in ${PDATA}/*.bam;
        do
        FILE_ID=`basename ${FILE} ".bam"`;

        SCRIPT=picard_${FILE_ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};	
        echo "java -Xmx4g -Djava.io.tmpdir=${TMP} \\
              -jar ${PCODE_PICARD}/SamToFastq.jar  \\
              INPUT=${FILE} \\
              FASTQ=${PDATA}/${FILE_ID}.fastq \\
              VERBOSITY=INFO" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/picard_${FILE_ID}.err -o ${PDATA}/picard_${FILE_ID}.out -v "BASH_ENV=~/.bashrc" -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
