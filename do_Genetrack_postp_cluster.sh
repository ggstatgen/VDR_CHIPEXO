#!/bin/bash
#postprocesses the output of do_Genetrack_cluster
#input: directory of gff files from Genetrack
#output: bed files with processed genetrack data
#for the kind of processing done, look into genetrack_postprocessing.pl

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input beds"
        exit
fi

PDATA=$1;
PCODE="/net/isi-scratch/giuseppe/scripts";

for FILE in ${PDATA}/*.gff;
        do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

	SCRIPT=${ID}.sh;

	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "perl ${PCODE}/genetrack_postprocessing.pl -input=${FILE} -d=35" >>${PDATA}/${SCRIPT};
        nice -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.bed -q medium_jobs.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
