#!/bin/bash
#script to perform sam validation with picard on the cluster
#Giuseppe Gallone 14/2/2013

#MODE=Mode	
#Mode of output Default value: VERBOSE. This option can be set to 'null' to clear the default value. 
#Possible values: {VERBOSE, SUMMARY}

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH_DATA> <EXTENSION>"
        echo "PATH_DATA - directory containing the input .sam files"
        echo "EXTENSION - [sam|bam]"
        exit
fi

PDATA=$1;
EXT=$2;

PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.117";
#human hg19 genome
#PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";
PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19_masked/hg19_masked.fa"
#temp
TMP="/tmp";


for FILE in ${PDATA}/*.${EXT};
        do
        FILE_ID=`basename ${FILE} ".${EXT}"`;

        SCRIPT=picard_validate_${FILE_ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        echo "java -Xmx4g -Djava.io.tmpdir=${TMP} \\
        -jar ${PCODE_PICARD}/ValidateSamFile.jar \\
        INPUT=${FILE} \\
        OUTPUT=${PDATA}/${FILE_ID}.samvalidation\\
        MODE=VERBOSE\\
        REFERENCE_SEQUENCE=${PFASTA} \\
        IGNORE_WARNINGS=true" >>${PDATA}/${SCRIPT};
      
        nice -5 qsub -e ${PDATA}/${FILE_ID}.stderr -o ${PDATA}/${FILE_ID}.stdout -v "BASH_ENV=~/.bashrc" -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
