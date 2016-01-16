#!/bin/bash

#24/9
#merge all bam files in a directory in order to do things on them afterwards (for example call phantompeak tools and see if the correlation produces meaningful peak pair distance)

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "PATH - path for the bam files (e.g. /home/me/files)"
        exit
fi

PDATA=$1;
TMP="/tmp";
PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.117";
SCRIPT=picard_bam_merge.sh;

for FILE in ${PDATA}/*.bam;
        do
        OPT="INPUT=";
        SPACE="  ";
        INPUTSTRING+=${OPT}${FILE}${SPACE};
        #need to build list of files
        #-I file1 -I file2 -I file3
done

echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};
echo "java -Xmx4g -Djava.io.tmpdir=${TMP} -jar ${PCODE_PICARD}/MergeSamFiles.jar \\
 ${INPUTSTRING} \\
OUTPUT=${PDATA}/picard_merged.sam \\
USE_THREADING=true" >>${PDATA}/${SCRIPT};
nice -5 qsub -e ${PDATA}/picard_bam_merge.err -o ${PDATA}/picard_bam_merge.out -q fgu217.q ${PDATA}/${SCRIPT};
#rm ${PDATA}/${SCRIPT};
