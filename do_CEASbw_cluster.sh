#!/bin/bash

#this is like CEAS however it accepts bigwig files, which you can easily obtain from MACS2 converting its bedgraph into bigwig
#downloaded from here
# https://bitbucket.org/cistrome/cistrome-applications-harvard/src/eda44af2697889584a18135c779cb3210ead2704/published-packages/CEAS?at=default


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
#  -l LENGTH_FILE, --length=LENGTH_FILE
#                        file contains lenth information of every chroms
#  --pf-res=PF_RES       Wig profiling resolution, DEFAULT: 50bp. WARNING:
#                        Value smaller than the wig interval (resolution) may
#                        cause aliasing error.





#This script assumes .wig and .bed have the same basename

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <TRACK_PATH> <WIG_RES>"
        echo "TRACK PATH - Where the .bed and .bw files are"
	echo "WIG_RES - wig sampling resolution. Default 50"
        exit
fi

PDATA=$1;
PRES=$2;

PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";
#Got this from the CEAS site: http://liulab.dfci.harvard.edu/CEAS/download.html
PGDB="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/CEAS/hg19.refGene";
#PGDB="/net/isi-scratch/giuseppe/indexes/Mmus/mm9/mm9.refGene"
PCHRS="/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes";
#PCHRS="/net/isi-scratch/giuseppe/indexes/chrominfo/mm9.chrom.sizes"

for FILE in ${PDATA}/*.bw;
	#do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
        do ID=`basename ${FILE} ".bw"`;

        SCRIPT=CEAS_${ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate' >>${PDATA}/${SCRIPT};
	#change required options here:
        echo "${PCODE}/ceasBW --bg -g ${PGDB} -b ${PDATA}/${ID}.bed -w ${FILE} --pf-res=${PRES} --name=${PDATA}/${ID} --length=${PCHRS}" >>${PDATA}/${SCRIPT};
 
        nice -5 qsub -e ${PDATA}/CEAS_${ID}.err -o ${PDATA}/CEAS_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
