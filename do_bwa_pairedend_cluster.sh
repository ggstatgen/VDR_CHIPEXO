#!/bin/bash
#script to run baw  with custom options on the cluster

#command: bwa aln -t8 -f [output.sai] [bwa_genome_index] [reads.fq]
#seguito da
#bwa samse [bwa_genome_index] [output.sai] [reads.fq] -f [output.sam]
#bwa sampe genome.index.fa s_4_1.sai s_4_2.sai s_4_1_sequence.fastq s_4_2_sequence.fastq > s_4.sam
#bwa sampe [options] <prefix> <in1.sai> <in2.sai> <in1.fq> <in2.fq>


#bwa aln options
#-e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]
#-l INT    seed length [32]
#-k INT    maximum differences in the seed [2]
#-I        the input is in the Illumina 1.3+ FASTQ-like format

#bwa sampe -f out.sam -r "@RQ\tID:<ID>\tLB:<LIBRARY_NAME>\tSM:<SAMPLE_NAME>\tPL:ILLUMINA" hg19 input1.sai input2.sai input1.fq input2.fq


PCODE="/net/isi-scratch/giuseppe/tools";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";
PINDEX="/net/isi-mirror/ucsc/mm10/chromosomes/mm10.fa";
PBWA="/home/giuseppe/local/bin" #this refers to a symbolic link to the latest bwa. If you don't use this, it will try to use the cgat bwa
#PINDEX="/net/isi-scratch/giuseppe/indexes/bwa/mm9.fa";
#qui ci vuole il mmus 9, crea un argomento linea di comando


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the fastq.gz files"
        exit
fi

PDATA=$1;

for FILE in ${PDATA}/*.fastq.gz;
	do
        PAIR=`echo ${FILE} | egrep -o "_2_"`;	
	if [ "${PAIR}" = '_2_' ]; 
	then
		continue;
	fi
 
        #BASEFILE=`basename ${FILE} "_1_sequence.fastq"`; 
	BASEFILE=`basename ${FILE} "_1_.fastq.gz"`;
        SCRIPT=${BASEFILE}.sh; #I should have 4 scripts only for this dataset
 
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
        echo "${PBWA}/bwa aln -q10 -t10 -f ${PDATA}/${BASEFILE}_1.sai ${PINDEX} ${PDATA}/${BASEFILE}_1_.fastq.gz" >>${PDATA}/${SCRIPT};
        echo "${PBWA}/bwa aln -q10 -t10 -f ${PDATA}/${BASEFILE}_2.sai ${PINDEX} ${PDATA}/${BASEFILE}_2_.fastq.gz" >>${PDATA}/${SCRIPT};
        echo "${PBWA}/bwa sampe ${PINDEX} ${PDATA}/${BASEFILE}_1.sai   ${PDATA}/${BASEFILE}_2.sai ${PDATA}/${BASEFILE}_1_.fastq.gz ${PDATA}/${BASEFILE}_2_.fastq.gz -f ${PDATA}/${BASEFILE}.sam" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "${PBWA}/samtools view -bS ${PDATA}/${BASEFILE}.sam > ${PDATA}/${BASEFILE}.bam" >>${PDATA}/${SCRIPT};
        #the stampy doc says stampy won't need sorted bams from bwa
        #comment this again
        echo "${PBWA}/samtools sort ${PDATA}${ID}.bam ${PDATA}${ID}.sorted" >>${PDATA}/${SCRIPT};
        echo "${PBWA}/samtools index ${PDATA}${ID}.sorted.bam" >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/${BASEFILE}.err -o ${PDATA}/${BASEFILE}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT}; 
	#find . -empty -type f -print0 | xargs -0 echo rm;
done
