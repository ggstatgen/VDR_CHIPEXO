#!/bin/bash

#modified from bgd2bw by Tao Liu
#needed to get bigwigs from the MACS2 output bedgraph because I need to feed them to ceasBW


if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <EXT> <CHROM_INFO>"
        echo "PATH - bedgraph directory"
	echo "EXT - file extension (eg. bdg)"
        echo "CHROM_INFO - file of chromosome sizes"
        exit
fi

PDATA=$1
EXT=$2
PCHROM=$3
PCODE="/net/isi-scratch/giuseppe/tools";
PBEDTOOLS="/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin";


#POUT=${PDATA/bdg/bw};
#echo ${POUT}
#mkdir ${POUT};

for FILE in ${PDATA}/*.${EXT};
        do BASEFILE=`basename ${FILE} ".${EXT}"`;
        ID=`basename ${FILE} ".${EXT}"`;

        SCRIPT=bdg2bw_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo 'source activate' >>${PDATA}/${SCRIPT};
	echo "${PBEDTOOLS}/slopBed -i ${FILE} -g ${PCHROM} -b 0 | ${PCODE}/UCSC_tools/bedClip stdin ${PCHROM} ${FILE}.clip" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        #echo "${PCODE}/UCSC_tools/bedGraphToBigWig ${FILE}.clip ${PCHROM} ${FILE/bdg/bw}" >>${PDATA}/${SCRIPT};
        echo "${PCODE}/UCSC_tools/bedGraphToBigWig ${FILE}.clip ${PCHROM} ${PDATA}/${BASEFILE}.bw" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "rm -f ${FILE}.clip" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/bdg2bw_${ID}.err -o ${PDATA}/bdg2bw_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done






# check commands: slopBed, bedGraphToBigWig and bedClip
#which slopBed &>/dev/null || which bedtools &>/dev/null || { echo "slopBed/bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
#which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
#which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
# end of checking

#if [ $# -lt 2 ];then
#    echo "Need 2 parameters! <bedgraph> <chrom info>"
#    exit
#fi
#
#F=$1
#G=$2
#
#${CODE}/bedtools-2.17.0/bin/slopBed -i ${F} -g ${G} -b 0 | ${CODE}/UCSC_tools/bedClip stdin ${G} ${F}.clip
#
#${CODE}/UCSC_tools/bedGraphToBigWig ${F}.clip ${G} ${F/bdg/bw}
#
#rm -f ${F}.clip
