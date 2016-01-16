#!/bin/bash
#script to perform sam->bam conversion and bam ordering using Picard, via FGU cluster
#Giuseppe Gallone 14/2/2013

############
#TASKS
############
#1. PICARD - sam -> bam conversion
#2. PICARD - mark PCR duplicates
#3. PICARD - fix Read Group information/generate index

#######################
#command line arguments
#######################
#1 data path (no final backslash)

#####
#INFO
#####
#The script uses the sam files' base filename as the seed name for all of the pipeline outputs
#INPUT: collection of SAM files obtained from an aligner of choice. BWA and BWA+STAMPY have been tested. All the .sam should be in the same directory
#all output files are currently stored in the INPUT directory


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input .sam files"
	exit
fi

PDATA=$1;

#------------------------
#hard coded paths - code
#------------------------
PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.84/picard-tools-1.84";
#------------------------
#hard coded paths - data
#------------------------
#human hg19 genome
PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19";
#temp
TMP="/tmp";

#------------------
#READ GROUP STRINGS
#-------------------
# !!! CHANGE HERE depending on dataset
#
RGID=""; #get this  dynamically from input file 
RGPL="ILLUMINA";
RGLB=""; #same as RGID for now
RGPU="STAMPY";  #change as needed
RGSM=""; #same as RGID for now
#
# !!! CHANGE HERE IF NEEDED
#
for FILE in ${PDATA}/*.sam;
        do
        FILE_ID=`basename ${FILE} ".sam"`;

        SCRIPT=sample_${FILE_ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        #----------------
	#Picard SAM->BAM
        #----------------
	echo 'echo "*START*: Picard SAM to BAM";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx4g -Djava.io.tmpdir=${TMP} \\
              -jar ${PCODE_PICARD}/SortSam.jar \\
              SORT_ORDER=coordinate \\
              INPUT=${FILE} \\
              OUTPUT=${PDATA}/${FILE_ID}.bam \\
              VALIDATION_STRINGENCY=LENIENT" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: Picard SAM to BAM";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

        #-----------------------
        #Picard Mark Duplicates
        #-----------------------
	echo 'echo "*START*: Picard Mark Duplicates";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx4g -Djava.io.tmpdir=${TMP}\\
              -jar ${PCODE_PICARD}/MarkDuplicates.jar\\
              INPUT=${PDATA}/${FILE_ID}.bam\\
              OUTPUT=${PDATA}/${FILE_ID}_marked.bam\\
              REMOVE_DUPLICATES=true\\
              METRICS_FILE=${PDATA}/metrics\\
              VALIDATION_STRINGENCY=LENIENT" >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: Picard Mark Duplicates";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	
	#----------------------
	#Picard Fix Read Group
	#----------------------
        echo 'echo "*START*: Picard Fix Read Group";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx4g -Djava.io.tmpdir=${TMP}\\
              -jar ${PCODE_PICARD}/AddOrReplaceReadGroups.jar\\
              INPUT=${PDATA}/${FILE_ID}_marked.bam\\
              OUTPUT=${PDATA}/${FILE_ID}_marked_rg.bam\\
              RGID=${FILE_ID}\\
              RGLB=${FILE_ID}\\
              RGPL=${RGPL}\\
              RGPU=${RGPU}\\
              RGSM=${FILE_ID}\\
              CREATE_INDEX=true" >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: Picard Fix Read Group";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/${FILE_ID}.stderr -o ${PDATA}/${FILE_ID}.stdout -v "BASH_ENV=~/.bashrc" -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
