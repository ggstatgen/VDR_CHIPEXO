#!/bin/bash

#This runs CNVnator on a bunch of whole genome bam sequences
#This is to produce Read Depth data as a prerequisite step for running the 
#alleleseq pipeline
#http://info.gersteinlab.org/AlleleSeq
#####
#(d) CNV file; a set of normalized read depth values for all snp locations from a separate genomic sequencing experiment. This reports the read depth at that snp compared to overall coverage. This is used to filter out locations with very low or high coverage, which would tend to indicate copy number variation. The file should be in this format (with header):
#
#chrm    snppos  rd
#1       52066   0.902113
#1       695745  0.909802
#1       742429  0.976435
#
#This can be generated from CNVnator (Abyzov et al. 2011). 

#This is an interesting blog post
#http://avrilomics.blogspot.co.uk/2013/01/using-cnvnator-to-find-copy-number.html

#CNVnator ASSUMES all the chromosomes fasta files ARE IN THE SAME DIRECTORY, OR A LINK IS PROVIDED

#This pipeline does the following:
#1
#>>>EXTRACTING READ MAPPING FROM BAM/SAM FILES
#$ ./cnvnator [-genome name] -root out.root [-chrom name1 ...] -tree [file1.bam ...]
#this creates a root object
#eg
#/net/isi-scratch/giuseppe/tools/CNVnator_v0.2.7/src/cnvnator -root NA19189.root -tree NA19189.mapped.ILLUMINA.bwa.YRI.low_coverage.20120522.bam

#2
#>>>GENERATING HISTOGRAM
#$ ./cnvnator [-genome name] -root file.root [-chrom name1 ...] -his bin_size [-d dir]
#eg
#/net/isi-scratch/giuseppe/tools/CNVnator_v0.2.7/src/cnvnator -genome /net/isi-scratch/giuseppe/indexes/Hsap/hs37d5.fa -root NA19189.root -his 100

#3
#>>>CALCULATING STATISTICS
#$ ./cnvnator -root file.root [-chrom name1 ...] -stat bin_size
#eg
#/net/isi-scratch/giuseppe/tools/CNVnator_v0.2.7/src/cnvnator -genome /net/isi-scratch/giuseppe/indexes/Hsap/hs37d5.fa -root NA19189.root -stat 100

#4
#>>>RD SIGNAL PARTITIONING
#$ ./cnvnator -root file.root [-chrom name1 ...] -partition bin_size [-ngc]
#eg
#/net/isi-scratch/giuseppe/tools/CNVnator_v0.2.7/src/cnvnator -genome /net/isi-scratch/giuseppe/indexes/Hsap/hs37d5.fa -root NA19189.root -partition 100

#>>>CNV CALLING
#$ ./cnvnator -root file.root [-chrom name1 ...] -call bin_size [-ngc]
#eg
#/net/isi-scratch/giuseppe/tools/CNVnator_v0.2.7/src/cnvnator -genome /net/isi-scratch/giuseppe/indexes/Hsap/hs37d5.fa -root NA19189.root -call 100 > cnvnator_NA19189_bin100_calls.tsv

#inputs: bin size
#RECOMMENDED:
#They recommend to use ~100-bp bins for 20-30x coverage, ~500-bp bins for 4-6x coverage, and ~30-bp bins for 100x coverage
#The bin size used shouldn't be shorter than the read length in your data.

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <WINDOW_SIZE>"
        echo "<PATH> - absolute path for genome sequencing bam files and for .fa chromosomes"
	echo "<BIN> - bin size - will depend on coverage, see scripts (eg: 100, 500,..)"
        exit
fi

PINPUT=$1;
BIN=$2;

#PCODE="/net/isi-scratch/giuseppe/tools/CNVnator_v0.2.7/src"; #path to cnvnator
PCODE="/net/isi-scratch/giuseppe/tools/CNVnator_v0.3/src"

#PGENOME="/net/isi-scratch/giuseppe/indexes/Hsap/hs37d5.fa";
PGENOME="/net/isi-scratch/giuseppe/indexes/Hsap/g1k_v37/g1k_v37_1_XYM.fa";


PDATA="/net/isi-scratch/giuseppe/VDR/ALLELESEQ/LC_16/CNV/d_root";

for FILE in ${PINPUT}/*.bam;
        do 
	ID=`basename ${FILE} ".bam"`;
	#eg  NA06986.mapped.ILLUMINA.bwa.CEU.low_coverage.20130415.bam
	#ID=`echo ${FILE} | egrep -o "NA[0-9]*"`;

        SCRIPT=CNVnator_${ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate'  >>${PDATA}/${SCRIPT};

	#export paths, just to be sure
	echo 'export ROOTSYS=/net/isi-scratch/giuseppe/tools/root' >>${PDATA}/${SCRIPT};
	echo 'export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib' >>${PDATA}/${SCRIPT};
	#change required options here:
	echo "cd ${PDATA}" >>${PDATA}/${SCRIPT};       
	#---------
	#1 - EXTRACTING READ MAPPING FROM BAM/SAM FILES
	#--------
	echo "${PCODE}/cnvnator -genome ${PGENOME} -root ${PDATA}/${ID}.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y  -tree ${FILE}" >>${PDATA}/${SCRIPT};
        #---------
        #2 - GENERATING HISTOGRAM
        #--------
	echo "${PCODE}/cnvnator -genome ${PGENOME} -root ${PDATA}/${ID}.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -his ${BIN}" >>${PDATA}/${SCRIPT};
        #---------
        #3 - CALCULATING STATISTICS
        #--------
        echo "${PCODE}/cnvnator -root ${PDATA}/${ID}.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y  -stat ${BIN}" >>${PDATA}/${SCRIPT};
        #---------
        #3 - RD SIGNAL PARTITIONING
        #--------
        echo "${PCODE}/cnvnator -root ${PDATA}/${ID}.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -partition ${BIN}" >>${PDATA}/${SCRIPT};
        #---------
        #4 - CNV CALLING
        #--------
        echo "${PCODE}/cnvnator -root ${PDATA}/${ID}.root  -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -call ${BIN} > ${PDATA}/cnvnator_${ID}_bin${BIN}_calls.tsv" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/CNVnator_${ID}.err -o ${PDATA}/CNVnator_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
