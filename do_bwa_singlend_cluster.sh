#!/bin/bash
#script to run baw  with custom options on the cluster

#command: bwa aln -t8 -f [output.sai] [bwa_genome_index] [reads.fq]
#seguito da
#bwa samse [bwa_genome_index] [output.sai] [reads.fq] -f [output.sam]
#bwa sampe genome.index.fa s_4_1.sai s_4_2.sai s_4_1_sequence.fastq s_4_2_sequence.fastq > s_4.sam
#bwa sampe [options] <prefix> <in1.sai> <in2.sai> <in1.fq> <in2.fq>


PCODE="/net/isi-scratch/giuseppe/tools";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Mmus/mm10/mm10.fa";
PINDEX="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Hsap/hg19_masked/hg19_masked.fa";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Hsap/hg19_nosex/hg19_nosex";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Hsap/g1k_v37/human_g1k_v37.fasta"
PBWA="/home/giuseppe/local/bin" #this refers to a symbolic link to the latest bwa. If you don't use this, it will try to use the cgat bwa
#PINDEX="/net/isi-scratch/giuseppe/indexes/bwa/mm9.fa";
#qui ci vuole il mmus 9, crea un argomento linea di comando
PSAMTOOLS="/home/giuseppe/local/bin/";


if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <PATH_DATA> <PATH_OUT> <LANE_ID> <QUEUE>"
        echo "PATH_DATA - directory containing the fastq files"
	echo "PATH_OUT - directory where to put the output bam files"
        echo "LANE_ID - lane name string to use for the Read Group PU field"
	echo "QUEUE - queue to send the job list to [newnodes|fgu217|medium_jobs]"
        exit
fi

PDATA=$1;
POUT=$2;
LANE=$3;
QUEUE=$4;

#------------------
#READ GROUP STRINGS
#-------------------
# !!! CHANGE HERE depending on dataset
#
RGPL="ILLUMINA";  #platform: gatk wants only one of ILLUMINA,SLX,SOLEXA,SOLID,454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN
RGLB="FOXN1"; #library 
RGPU=${LANE};  #platform unit (lane info) 
RGPG="BWA"; #Programs used for processing the read group.
RGCN="Peconic";
#
# !!! CHANGE HERE AS NEEDED
#

for FILE in ${PDATA}/*.fastq;
	do 
        BASEFILE=`basename ${FILE} ".fastq"`; 
        #FILE_ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;
	FILE_ID=`basename ${FILE} ".fastq"`;

	#sample dependent variables
	SCRIPT=bwa_se_${FILE_ID}.sh;
	RGSM=${FILE_ID};
	#ID: Should be a combination of sample (SM), library (LB), lane (PU)
	RGID=${RGLB}_${RGSM}_${RGPU};

        echo '#!/bin/bash' >>${POUT}/${SCRIPT};
        echo '' >>${POUT}/${SCRIPT};

	#change required options here:
        echo "${PBWA}/bwa aln -q10 -t15 -f ${POUT}/${BASEFILE}.sai ${PINDEX} ${FILE}" >>${POUT}/${SCRIPT};
        echo "${PBWA}/bwa samse -r \"@RG\tID:${RGID}\tLB:${RGLB}\tSM:${RGSM}\tPL:${RGPL}\tPG:${RGPG}\tCN:${RGCN}\tPU:${RGPU}\"  ${PINDEX} ${POUT}/${BASEFILE}.sai ${FILE} -f ${POUT}/${BASEFILE}.sam" >>${POUT}/${SCRIPT};
        echo '' >>${POUT}/${SCRIPT};
	echo "${PSAMTOOLS}/samtools view -bS ${POUT}/${BASEFILE}.sam > ${POUT}/${BASEFILE}.bam" >>${POUT}/${SCRIPT};
        
       	nice -5 qsub -e ${POUT}/bwa_se_${BASEFILE}.err -o ${POUT}/bwa_se_${BASEFILE}.out -q ${QUEUE}.q ${POUT}/${SCRIPT};
        rm ${POUT}/${SCRIPT}; 
	#find . -empty -type f -print0 | xargs -0 echo rm;
done
