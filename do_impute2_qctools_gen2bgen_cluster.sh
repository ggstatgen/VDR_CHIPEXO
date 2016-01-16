#!/bin/bash

#this converts text .gen files produced by IMPUTE2 in binary .bgen files
#these are smaller and I use them as input for SNPTEST

#usage
#http://www.well.ox.ac.uk/~gav/qctool/#tutorial
#qctool -g example_#.gen -og example_#.bgen


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input .gen files"
        exit
fi

PDATA=$1;
PCODE="/net/isi-scratch/giuseppe/tools/qctool_v1.3-linux-x86_64";

for FILE in ${PDATA}/*.gen;
        do ID=`basename ${FILE} ".gen"`;
	
	SCRIPT=qctools_gen2bgen_${ID}.sh;

	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "${PCODE}/qctool -g ${FILE} -og ${PDATA}/${ID}.bgen;" >>${PDATA}/${SCRIPT};
	nice -5 qsub -e ${PDATA}/qctools_gen2bgen_${ID}.err -o ${PDATA}/qctools_gen2bgen_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
