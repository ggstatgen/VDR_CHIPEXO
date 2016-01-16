#!/bin/bash
#script to run baw  with custom options on the cluster

#command: bwa aln -t8 -f [output.sai] [bwa_genome_index] [reads.fq]
#seguito da
#bwa samse [bwa_genome_index] [output.sai] [reads.fq] -f [output.sam]
#bwa sampe genome.index.fa s_4_1.sai s_4_2.sai s_4_1_sequence.fastq s_4_2_sequence.fastq > s_4.sam
#bwa sampe [options] <prefix> <in1.sai> <in2.sai> <in1.fq> <in2.fq>


PCODE="/net/isi-scratch/giuseppe/tools";
PINDEX="/net/isi-scratch/giuseppe/indexes/Mmus/mm10/mm10.fa";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Hsap/hg19_masked/hg19_masked.fa";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Hsap/hg19_nosex/hg19_nosex";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Hsap/g1k_v37/human_g1k_v37.fasta"
PBWA="/home/giuseppe/local/bin" #this refers to a symbolic link to the latest bwa. If you don't use this, it will try to use the cgat bwa
#PINDEX="/net/isi-scratch/giuseppe/indexes/bwa/mm9.fa";
#qui ci vuole il mmus 9, crea un argomento linea di comando

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the fastq files"
        exit
fi

PDATA=$1;
POUT=${PDATA}/d_BWA;
mkdir ${POUT};

for FILE in ${PDATA}/*.fastq;
	do 
        ID=`basename ${FILE} ".fastq"`; 

	#sample dependent variables
	SCRIPT=bwa_se_${ID}.sh;
        echo '#!/bin/bash' >>${POUT}/${SCRIPT};
        echo '' >>${POUT}/${SCRIPT};

	#change required options here:
        echo "${PBWA}/bwa aln -q10 -t15 -f ${POUT}/${ID}.sai ${PINDEX} ${FILE}" >>${POUT}/${SCRIPT};
        echo "${PBWA}/bwa samse ${PINDEX} ${POUT}/${ID}.sai ${FILE} -f ${POUT}/${ID}.sam" >>${POUT}/${SCRIPT};
        echo '' >>${POUT}/${SCRIPT};
	echo "${PCODE}/samtools/samtools view -bS ${POUT}/${ID}.sam > ${POUT}/${ID}.bam" >>${POUT}/${SCRIPT};
        
       	nice -5 qsub -e ${POUT}/bwa_se_${ID}.err -o ${POUT}/bwa_se_${ID}.out -q newnodes.q ${POUT}/${SCRIPT};
        rm ${POUT}/${SCRIPT}; 
done
