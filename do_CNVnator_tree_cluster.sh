#!/bin/bash

#This runs CNVnator on a bunch of fully sequenced genomes
#This is to produce Read Depth data as a prerequisite step for running the 
#alleleseq pipeline
#http://info.gersteinlab.org/AlleleSeq
#####

#This version only gets the basic root file and will be the input for the 
#scripts provided by the guys at the gerstein lab
#see mail 3/2/14 Jieming Chen

#>>>EXTRACTING READ MAPPING FROM BAM/SAM FILES
#
#$ ./cnvnator [-genome name] -root out.root [-chrom name1 ...] -tree [file1.bam ...]
#
#out.root  -- output ROOT file. See ROOT package documentation.
#chr_name1 -- chromosome name.
#file.bam  -- bam files.




if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "<PATH> - absolute path for genome sequencing bam files and for .fa chromosomes"
        exit
fi

PDATA=$1;

#PCODE="/net/isi-scratch/giuseppe/tools/CNVnator_v0.2.7/src"; #path to cnvnator
PCODE="/net/isi-scratch/giuseppe/tools/CNVnator_v0.3/src";

#PGENOME="/net/isi-scratch/giuseppe/indexes/Hsap/hs37d5.fa";
PGENOME="/net/isi-scratch/giuseppe/indexes/Hsap/g1k_v37/g1k_v37_1_XYM.fa";
#PGENOME="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";

for FILE in ${PDATA}/NA*.bam;
        do 
	ID=`basename ${FILE} ".bam"`;
	#ID=`echo ${FILE} | egrep -o "NA[0-9]*.mapped"`;

        SCRIPT=CNVnator_${ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate'  >>${PDATA}/${SCRIPT};
	#export paths, just to be sure
#	echo 'ROOTSYS=/net/isi-scratch/giuseppe/tools/root' >>${PDATA}/${SCRIPT};
	echo 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/net/isi-scratch/giuseppe/tools/root/lib' >>${PDATA}/${SCRIPT};
	#change required options here:
	echo "cd ${PDATA}" >>${PDATA}/${SCRIPT};       
	#---------
	#1 - EXTRACTING READ MAPPING FROM BAM/SAM FILES
	#--------
	#echo "${PCODE}/cnvnator -root ${ID}.root -tree ${FILE}" >>${PDATA}/${SCRIPT};
	echo "${PCODE}/cnvnator -genome ${PGENOME} -root ${ID}.root  -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -tree ${FILE}" >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/CNVnator_${ID}.err -o ${PDATA}/CNVnator_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        #nice -5 qsub -v "BASH_ENV=~/.bashrc" -e ${PDATA}/CNVnator_${ID}.err -o ${PDATA}/CNVnator_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
