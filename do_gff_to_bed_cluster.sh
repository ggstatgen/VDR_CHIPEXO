#!/bin/bash

#uses bedops
#Usage:
#  /net/isi-scratch/giuseppe/tools/bedops/bin/gff2bed [ --help ] [ --do-not-sort | --max-mem <value> ] < foo.gff
#
#Options:
#  --help              Print this help message and exit
#  --do-not-sort       Do not sort converted data with BEDOPS sort-bed
#  --max-mem <value>   Sets aside <value> memory for sorting BED output. For example,
#                     <value> can be 8G, 8000M or 8000000000 to specify 8 GB of memory
#                     (default: 2G).
#
#About:
#  This script converts 1-based, closed [a,b] GFF3 data from standard
#  input into 0-based, half-open [a-1,b) six-column extended BED, sorted and
#  sent to standard output.
#
#  The GFF3 specification (http://www.sequenceontology.org/gff3.shtml)
#  contains columns that do not map directly to common or UCSC BED columns.
#  Therefore, we add the following columns to preserve the ability to
#  seamlessly convert back to GFF3 after performing operations with
#  bedops, bedmap, or other BEDOPS or BED-processing tools.
#
#  - The 'source' GFF column data maps to the 7th BED column
#  - The 'type' data maps to the 8th BED column
#  - The 'phase' data maps to the 9th BED column
#  - The 'attributes' data maps to the 10th BED column
#
#  We make the following assumptions about the GFF3 input data:
#
#  - The 'seqid' GFF column data maps to the chromosome label (1st BED column)
#  - The 'ID' attribute in the 'attributes' GFF column (if present) maps to
#    the element ID (4th BED column)
#  - The 'score' and 'strand' GFF columns (if present) are mapped to the
#    5th and 6th BED columns, respectively
#
#  If we encounter zero-length insertion elements (which are defined
#  where the start and stop GFF column data values are equivalent), the
#  start coordinate is decremented to convert to 0-based, half-open indexing,
#  and a 'zero_length_insertion' attribute is added to the 'attributes' GFF
#  column data.


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "You must specify the data path for the gff files (e.g. /home/me/files)"
        exit
fi

PDATA=$1;
PCODE="/net/isi-scratch/giuseppe/tools/bedops/bin";

for FILE in ${PDATA}/*.gff;
	do ID=`basename ${FILE} ".gff"`;

        SCRIPT=gff2bed_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	
	echo "source activate" >> ${PDATA}/${SCRIPT};
	echo "${PCODE}/gff2bed --do-not-sort < ${FILE} | cut -f 1,2,3,4 | sort -k1,1V -k2,2g > ${PDATA}/${ID}.clean.bed" >> ${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/gff2bed_${ID}.err -o ${PDATA}/gff2bed_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
