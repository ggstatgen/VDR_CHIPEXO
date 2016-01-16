#!/bin/bash
#uses the genetrack peakcaller to find peaks in .bed files obtained by converting bam files using bedtools

#Usage: genetrack.py [options] input_paths
#
#input_paths may be:
#- a file to run on
#- "-" to run on standard input
#
#example usages:
#python genetrack.py -s 10 /path/to/a/file.txt
#python genetrack.py -s 5 -e 50 -
#
#Options:
#  -h, --help     show this help message and exit
#  -s SIGMA       Sigma to use when smoothing reads to call peaks. Default 5
#  -e EXCLUSION   Exclusion zone around each peak that prevents others from
#                 being called. Default 20.
#  -u UP_WIDTH    Upstream width of called peaks. Default uses half exclusion
#                 zone.
#  -d DOWN_WIDTH  Downstream width of called peaks. Default uses half exclusion
#                 zone.
#  -F FILTER      Absolute read filter; outputs only peaks with larger peak
#                 height. Default 3.
#  -c CHROMOSOME  Chromosome (ex chr11) to limit to. Default process all.
#  -k CHUNK_SIZE  Size, in millions of base pairs, to chunk each chromosome
#                 into when processing. Each 1 million size uses approximately
#                 20MB of memory. Default 10.
#  -o FORMAT      Output format for called peaks. Valid formats are gff and
#                 txt. Default gff.
#  -b             Output bed graph tracks.
#  -v             Verbose mode: displays debug messages

#NOTE 2012/12/10 I modified genetrack.py to rename output bedgraphs from "forward/reverse.bedgraph" to "inputnamefile_forward/reverse.bedgraph"
#TODO: modify the script to output in the current dir


#NOTE The peak calling section of the peconics documentations says they do two runs:
#fine grain -s 5 -e 10
#coarse grain -s 10 -e 40

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input beds"
        exit
fi

PDATA=$1;
#PCODE="/net/isi-scratch/giuseppe/tools/chipexo-master/genetrack";
PCODE="/net/isi-scratch/giuseppe/tools/chipexo-master_modded_GG/genetrack"; #use this if you want the bedtrack files

for FILE in ${PDATA}/*.bed;
        do 
	#ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
	ID=`basename ${FILE} ".bed"`;

	SCRIPT=genetrack_${ID}.sh;

	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate' >>${PDATA}/${SCRIPT};
	echo "python ${PCODE}/genetrack.py -v -s 5 -e 10 -b ${FILE}" >>${PDATA}/${SCRIPT};
	#echo "python ${PCODE}/genetrack.py -v -s 5 -e 10 -F 6 ${FILE}" >>${PDATA}/${SCRIPT};
        nice -5 qsub -e ${PDATA}/genetrack_${ID}.err -o ${PDATA}/genetrack_${ID}.gff -q newnodes.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
