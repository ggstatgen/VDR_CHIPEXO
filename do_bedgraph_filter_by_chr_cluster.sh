#!/bin/bash

#2/6/2014
#Trying to visualise normalised signal tracks in igv or UCSC
#seems like both complain they're too big
#also converted bw from these tracks won't fit into my dropbox folder

#here I split them by chromosomes, based on the interesting peaks I want to visualise


if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <EXT> <SEARCH_STRING>"
        echo "<PATH> data path for the bedgraph (or bed file)"
	echo "<EXT> file extension"
	echo "<SEARCH_STRING> what to grep for (e.g. chr8)"
        exit
fi

PDATA=$1;
EXT=$2;
GREP=$3;

for FILE in ${PDATA}/*.${EXT};
        #do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
        do ID=`basename ${FILE} ".${EXT}"`;

        SCRIPT=sel_${GREP}_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "grep -P  \"^${GREP}\W\" ${FILE} > ${PDATA}/${GREP}_${ID}.${EXT}" >> ${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/sel_${GREP}_${ID}.err -o ${PDATA}/sel_${GREP}_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
