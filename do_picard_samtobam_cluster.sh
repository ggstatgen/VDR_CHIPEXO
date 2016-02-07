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
PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.117";
#------------------------
#hard coded paths - data
#------------------------
#human hg19 genome
#PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19";
#PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/g1k_v37/human_g1k_v37.fasta"
#temp
TMP="/tmp";

#------------------
#READ GROUP STRINGS
#-------------------
# !!! CHANGE HERE depending on dataset
#
RGID="MATTEA_CHIPEXO" #append to this dynamically from input file 
RGPL="ILLUMINA";  #gatk wants only one of ILLUMINA,SLX,SOLEXA,SOLID,454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN
RGLB="-";
RGPU="BWA_mm10";  #really don't know what to put in this one
#RGSM=""; #Giulio set this to "Andrew" and "Jean", I cannot get these from the file names, will use some variant of RGID taken from the file name
#
# !!! CHANGE HERE IF NEEDED
#


for FILE in ${PDATA}/*.sam;
        do
        ID=`basename ${FILE} ".sam"`;

        SCRIPT=picard_s2b_${ID}.sh;
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
              OUTPUT=${PDATA}/${ID}.bam \\
              VALIDATION_STRINGENCY=LENIENT" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: Picard SAM to BAM";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

        #-----------------------
        #Picard Mark Duplicates
        #-----------------------
	#removed: REMOVE_DUPLICATES=true\\
	echo 'echo "*START*: Picard Mark Duplicates";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx4g -Djava.io.tmpdir=${TMP}\\
              -jar ${PCODE_PICARD}/MarkDuplicates.jar\\
              INPUT=${PDATA}/${ID}.bam\\
              OUTPUT=${PDATA}/${ID}_marked.bam\\
              METRICS_FILE=${PDATA}/${ID}.metrics\\
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
              INPUT=${PDATA}/${ID}_marked.bam\\
              OUTPUT=${PDATA}/${ID}_final.bam\\
              RGID=${RGID}_${ID}\\
              RGLB=${RGLB}\\
              RGPL=${RGPL}\\
              RGPU=${RGPU}\\
              RGSM=sample_${ID}\\
              CREATE_INDEX=true" >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: Picard Fix Read Group";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/picard_s2b_${ID}.err -o ${PDATA}/picard_s2b_${ID}.out -v "BASH_ENV=~/.bashrc" -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
