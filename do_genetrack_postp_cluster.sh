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
PCODE="/net/isi-backup/giuseppe/scripts";

for FILE in ${PDATA}/*.gff;
        do 
	ID=`basename ${FILE} ".gff"`;
	SCRIPT=gt_pp_${ID}.sh;

	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "perl ${PCODE}/genetrack_postprocessing.pl -input=${FILE} -d=40" >>${PDATA}/${SCRIPT};
        nice -5 qsub -e ${PDATA}/gt_pp_${ID}.err -o ${PDATA}/${ID}.bed -q newnodes.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
