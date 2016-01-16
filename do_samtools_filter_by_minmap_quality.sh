#!/bin/bash

#following the samtools faq 
#http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_FAQ#I_want_to_get_.60unique.27_alignments_from_SAM.2FBAM.
#this script should remove all bwa/stampy mappings which map to non-unique location, or better which have 0 quality.

#samtools view -q INT   minimum mapping quality [0]

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <MIN_MAP_QUAL>"
        echo "You must specify the data path for the bam files (e.g. /home/me/files) and a number for the minimum acceptable mapping quality"
        exit
fi

PDATA=$1;
PCODE="/home/giuseppe/local/bin";

for FILE in ${PDATA}/*.bam;
	#change this depending on the data
	#do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
        do 
        #ID=`echo ${FILE} | egrep -o "38380_[0-9]*"`;        
        
        ID=`basename ${FILE} ".bam"`;
        
        SCRIPT=st_minmap_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
        echo "${PCODE}/samtools view -bq $2 ${FILE} > ${PDATA}/${ID}-minmapq_$2.bam" >>${PDATA}/${SCRIPT};
        nice -5 qsub -v "BASH_ENV=~/.bashrc" -e ${PDATA}/st_minmap_${ID}.err -o ${PDATA}/st_minmap_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
