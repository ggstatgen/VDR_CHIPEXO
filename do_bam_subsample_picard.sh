#!/bin/bash

#script to use picard to subsample a bam in order to equalise sizes
#This is just a basic script to see if all works correctly
#a better approach would be to subsample many times and call peaks from x random samples, and only keep peaks shared by all samples

#picard DownsampleSam
#Randomly down-sample a SAM or BAM file to retain a random subset of the reads. 
#Mate-pairs are either both kept or both discarded. 
#Reads marked as not primary alignments are all discarded. 
#Each read is given a probability P of being retained - results with the exact same input 
#in the same order and with the same value for RANDOM_SEED will produce the same results.

#Option	Description
#INPUT=File	The input SAM or BAM file to downsample. Required.
#OUTPUT=File	The output, downsampled, SAM or BAM file to write. Required.
#RANDOM_SEED=Long	Random seed to use if reproducibilty is desired. Setting to null will cause multiple invocations to produce different results. Default value: 1. This option can be set to 'null' to clear the default value.
#PROBABILITY=Double	The probability of keeping any individual read, between 0 and 1. Default value: 1.0. This option can be set to 'null' to clear the default value. 


if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH_DATA> <PROBABILITY>"
        echo "PATH_DATA - directory containing the input .bam files"
	echo "PROBABILITY - The probability of keeping an individual read, from 0 to 1"
        exit
fi

PDATA=$1;
PROB=$2;

PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.102";
TMP="/tmp";

for FILE in ${PDATA}/*.bam;
        do
        ID=`basename ${FILE} ".bam"`;

        SCRIPT=picard_subsample_${ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx2g -Djava.io.tmpdir=${TMP} \\
              -jar ${PCODE_PICARD}/DownsampleSam.jar \\
              INPUT=${FILE} \\
              OUTPUT=${PDATA}/${ID}_downsampled.bam \\
              PROBABILITY=${PROB} \\
              VALIDATION_STRINGENCY=LENIENT" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	
	nice -5 qsub -e ${PDATA}/picard_subsample_${ID}.err -o ${PDATA}/picard_subsample_${ID}.out -v "BASH_ENV=~/.bashrc" -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
