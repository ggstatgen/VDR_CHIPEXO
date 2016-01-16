#!/bin/bash
#to index bgzipped vcf files using tabix

#Usage:   tabix <in.tab.bgz> [region1 [region2 [...]]]
#
#Options: -p STR     preset: gff, bed, sam, vcf, psltbl [gff]
#         -s INT     sequence name column [1]
#         -b INT     start column [4]
#         -e INT     end column; can be identical to '-b' [5]
#         -S INT     skip first INT lines [0]
#         -c CHAR    symbol for comment/meta lines [#]
#         -r FILE    replace the header with the content of FILE [null]
#         -B         region1 is a BED file (entire file will be read)
#         -0         zero-based coordinate
#         -h         print the header lines
#         -l         list chromosome names
#         -f         force to overwrite the index

if [ ! $# == 2 ]; then
	echo "Usage: `basename $0` <PATH> <EXTENSION>"
        echo "You must specify 1) data path (e.g. /home/me/files) and 2) file extension (e.g. vcf.gz)"
	exit
fi

PDATA=$1;
PCODE="/net/isi-scratch/giuseppe/tools/tabix-0.2.6";
EXT=$2;

for FILE in ${PDATA}/*.${EXT};
        #do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
	do ID=`basename ${FILE} ".${EXT}"`;	

	SCRIPT=${ID}.sh;
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "${PCODE}/tabix -f ${FILE}" >>${PDATA}/${SCRIPT};
	nice -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
