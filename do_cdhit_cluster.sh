#!/bin/bash

#CD-HIT clusters proteins into clusters that meet a user-defined similarity threshold, usually a sequence identity. Each cluster has one representative sequence. The input is a protein dataset in fasta format and the output are two files: a fasta file of representative sequences and a text file of list of clusters.

#CD-HIT-EST

#CD-HIT-EST clusters a nucleotide dataset into clusters that meet a user-defined similarity threshold, usually a sequence identity. The input is a DNA/RNA dataset in fasta format and the output are two files: a fasta file of representative sequences and a text file of list of clusters. Since eukaryotic genes usually have long introns, which cause long gaps, it is difficult to make full-length alignments for these genes. So, CD-HIT-EST is good for non-intron containing sequences like EST.

#Basic command:

#cd-hit-est -i est_human -o est_human95 -c 0.95 -n 8  
#Choose of word size:

#-n 8,9,10 for thresholds 0.90 ~ 1.0
#-n 7      for thresholds 0.88 ~ 0.9
#-n 6      for thresholds 0.85 ~ 0.88
#-n 5      for thresholds 0.80 ~ 0.85
#-n 4      for thresholds 0.75 ~ 0.8 


PCODE="/net/isi-scratch/giuseppe/tools";

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the fasta files"
        exit
fi

PDATA=$1;

for FILE in ${PDATA}/*.fasta;
	do 
        BASEFILE=`basename ${FILE} ".fasta"`; 
	SCRIPT=${BASEFILE}.sh;
 
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	echo "${PCODE}/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i ${FILE} -o ${PDATA}/${BASEFILE} -M 0 -T 15" >>${PDATA}/${SCRIPT};
       	nice -5 qsub -e ${PDATA}/${BASEFILE}.err -o ${PDATA}/${BASEFILE}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT}; 
done
