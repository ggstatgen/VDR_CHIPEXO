#!/bin/bash

#CEAS - Cis-regulatory Element Annotation System
#http://liulab.dfci.harvard.edu/CEAS/usermanual.html
#Usage: ceas < input files > [options]
#basic options:
#  -b BED, --bed=BED     BED file of ChIP regions.
#  -w WIG, --wig=WIG     WIG file for either wig profiling or genome background
#                        annotation. WARNING: --bg flag must be set for genome
#                        background re-annotation.
#  -e EBED, --ebed=EBED  BED file of extra regions of interest (eg, non-coding
#                        regions)
#  -g GDB, --gt=GDB      Gene annotation table (eg, a refGene table in sqlite3
#                        db format provided through the CEAS web,
#                        http://liulab.dfci.harvard.edu/CEAS/download.html).

#This script assumes .wig and .bed have the same basename

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <TRACK_PATH>"
        echo "TRACK PATH - Where the .bed and .wig files are"
        exit
fi

PDATA=$1;
PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";
#Got this from the CEAS site: http://liulab.dfci.harvard.edu/CEAS/download.html
PGDB="/net/isi-scratch/giuseppe/indexes/Hsap/CEAS/hg19.refGene";

for FILE in ${PDATA}/*.bed;
	#do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
        do ID=`basename ${FILE} ".bed"`;

        SCRIPT=${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
        #echo "${PCODE}/ceas -b ${FILE} -w ${PDATA}/${ID}.wig -g ${PGDB} --name=${PDATA}/${ID}" >>${PDATA}/${SCRIPT};
        echo "${PCODE}/ceas -b ${FILE} -g ${PGDB} --name=${PDATA}/${ID}" >>${PDATA}/${SCRIPT};

        
        nice -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
