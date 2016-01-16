#!/bin/bash
#input: sam file from bwa
#output: sam file with RG field for each read

#needed before using GATK - I want to tag each read with info, including lane info, THEN fuse the sam files for all lanes

if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH_DATA> <LANE_ID> <QUEUE>"
        echo "PATH_DATA - directory containing the fastq files"
        echo "LANE_ID - lane name string to use for the Read Group PU field"
        echo "QUEUE - queue to send the job list to [newnodes|fgu217|medium_jobs]"
        exit
fi

PDATA=$1;
LANE=$2;
QUEUE=$3;


#------------------------
#hard coded paths - code
#------------------------
PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.117";
TMP="/tmp"; #GATK gives "too many open files" exception here if you use too many threads

#------------------
#READ GROUP STRINGS
#-------------------
# !!! CHANGE HERE depending on dataset
#
RGPL="ILLUMINA";  #platform: gatk wants only one of ILLUMINA,SLX,SOLEXA,SOLID,454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN
RGLB="VDR_sc-1008_LCL"; #library 
RGPU=${LANE};  #platform unit (lane info) 
RGPG="BWA"; #Programs used for processing the read group. PICARD WONT RECOGNISE THIS
RGCN="Peconic";
#
# !!! CHANGE HERE AS NEEDED
#

for FILE in ${PDATA}/*.sam;
        do
        FILENAME=`basename ${FILE} ".sam"`;
	FILE_ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        SCRIPT=picard_fix_rg_${FILE_ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        RGSM=${FILE_ID};
        #ID: Should be a combination of sample (SM), library (LB), lane (PU)
        RGID=${RGLB}_${RGSM}_${RGPU};

        #----------------
        #Picard SAM->BAM
        #----------------
        #echo 'echo "*START*: Picard SAM to BAM";' >>${PDATA}/${SCRIPT};
        #echo '' >>${PDATA}/${SCRIPT};
        #echo "java -Xmx4g -Djava.io.tmpdir=${TMP} \\
        #      -jar ${PCODE_PICARD}/SortSam.jar \\
        #      SORT_ORDER=coordinate \\
        #      INPUT=${FILE} \\
        #      OUTPUT=${PDATA}/${FILENAME}.bam \\
        #      VALIDATION_STRINGENCY=LENIENT" >>${PDATA}/${SCRIPT};
        #echo '' >>${PDATA}/${SCRIPT};
        #echo 'echo "*DONE*: Picard SAM to BAM";' >>${PDATA}/${SCRIPT};
        #echo '' >>${PDATA}/${SCRIPT};
       
	#----------------------
        #Picard Fix Read Group
        #----------------------
        echo 'echo "*START*: Picard Fix Read Group";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx4g -Djava.io.tmpdir=${TMP}\\
              -jar ${PCODE_PICARD}/AddOrReplaceReadGroups.jar\\
              INPUT=${FILE}\\
              OUTPUT=${PDATA}/${FILENAME}_rg.sam\\
              RGSM=${RGSM}\\
              RGID=${RGID}\\
              RGLB=${RGLB}\\
              RGPL=${RGPL}\\
              RGPU=${RGPU}\\
              RGCN=${RGCN}\\
              CREATE_INDEX=false" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "*DONE*: Picard Fix Read Group";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/${FILE_ID}.err -o ${PDATA}/${FILE_ID}.out -v "BASH_ENV=~/.bashrc" -q ${QUEUE}.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
