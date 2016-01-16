#!/bin/bash
#script to run fastqc with custom options on the cluster

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH_DATA> <CHROMOSOMESIZE_FILE>"
        echo "PATH_DATA - directory containing the *compressed* wig files"
        echo "CHROMOSOMESIZE_FILE -  absolute path for the file of chromosome sizes (eg there are some in /net/isi-scratch/giuseppe/tools/UCSC_tools)"
        exit
fi

PDATA=$1;
PCHROM=$2;
#location of wigtobigwig executable, change if needed
PCODE="/net/isi-scratch/giuseppe/tools/UCSC_tools";

for FILE in ${PDATA}/*.wig;
	do 
        BASEFILE=`basename ${FILE} ".wig"`; 
        #ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
        SCRIPT=w_to_bw_${BASEFILE}.sh;
        
        #remove scripts if it already exists
        #if [ -f "filename" ]
	#then
        #rm ${PDATA}/${SCRIPT};
        #fi
 
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
        echo "${PCODE}/wigToBigWig -clip ${FILE} ${PCHROM} ${PDATA}/${BASEFILE}.bw;" >>${PDATA}/${SCRIPT};
        nice -19 qsub -e ${PDATA}/${BASEFILE}.err -o ${PDATA}/${BASEFILE}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT}; 

	#find . -empty -type f -print0 | xargs -0 echo rm;
done
