#!/bin/bash

#script to run the platypus variant caller using the cluster
#Running Platypus
#===================
#
#The easiest way to run Platypus is as follows:
#
#python Platypus.py callVariants --bamFiles=LIST_OF_BAMS --refFile=REF.fa --output=Calls.vcf
#
#
#You can see a list of all the possible input options by running the following comand:
#
#python Platypus.py callVariants --help
#
#However, in most cases the default parameter values should be fine, and you will only need to specify the --bamFiles
#and --refFile and --output arguments. By default, if you do not specify a region or regions or interest, Platypus will
#run through all the data in your BAM files. The --regions argument can be used to specify regions of interest.

#3. Running in Variant-Calling Mode
#==================================
#
#The standard way of running Platypus is to use it to detect variants in one or more BAM files. Variants are detected by
#comparing the BAM reads with a reference sequence. This can be done using the following command:
#
#python Platypus.py callVariants --bamFiles=DATA.bam --regions=chr20 --output=test.vcf --refFile=GENOME.fa
#
#where the input BAM files, and the genome reference must be indexed using samtools, or a program that produces compatible
#index files.

#--nCPU

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "You must specify the data path for the bam files (e.g. /home/me/files)"
        exit
fi

PDATA=$1;
#location of macs2 executable, change if needed
PCODE="/net/isi-scratch/giuseppe/tools/Platypus_0.5.2";
REFERENCE="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";

for FILE in ${PDATA}/*.bam;
	do ID=`basename ${FILE} ".bam"`;

        SCRIPT=platypus_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	echo "source activate" >> ${PDATA}/${SCRIPT};
	echo "python ${PCODE}/Platypus.py callVariants --nCPU 10 --bamFiles=${FILE} --refFile=${REFERENCE} --output=${PDATA}/platypus_${ID}_calls.vcf" >>${PDATA}/${SCRIPT};

        nice -5 qsub -pe dedicated 10 -e ${PDATA}/platypus_${ID}.err -o ${PDATA}/platypus_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};

        #rm ${PDATA}/${SCRIPT};  
done
