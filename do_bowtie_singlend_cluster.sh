#!/bin/bash
#script to run bowtie 2.19 with options taken from the GEM paper
#the idea is to use alignments for GEM which were obtained exactly as obtained for GEM

#command:
#bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]
#I will use bowtie -q --best --strata --m 1 -p 15  --chunkmbs 1024
#where
#  -q                 query input files are FASTQ .fq/.fastq (default)
#  --best             hits guaranteed best stratum; ties broken by quality
#  --strata           hits in sub-optimal strata aren't reported (requires --best)
#  -m <int>           suppress all alignments if > <int> exist (def: no limit)
#  -p/--threads <int> number of alignment threads to launch (default: 1)
#  --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 64)
#  -S                 sam output

PCODE="/net/isi-scratch/giuseppe/tools";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19_bowtie";
#PINDEX="/net/isi-scratch/giuseppe/indexes/Mmus/bwa/mm9_bowtie";
PINDEX="/net/isi-scratch/giuseppe/indexes/Mmus/mm10/mm10_bowtie"
#PINDEX="/net/isi-scratch/giuseppe/indexes/Hsap/g1k_v37/g1k_v37_1_XYM"

#PINDEX="/net/isi-scratch/giuseppe/indexes/"
#qui ci vuole il mmus 9, crea un argomento linea di comando


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the fastq files"
        exit
fi

PDATA=$1;

for FILE in ${PDATA}/*.fastq;
	do 

        ID=`basename ${FILE} ".fastq"`;
	#ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        BASEFILE=`basename ${FILE} ".fastq"`; 
	SCRIPT=bowtie_se_${BASEFILE}.sh;
 
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
	#this is for Alleleseq (same option as in PIPELINE.mk)
	#--best --strata -v 2 -m 1 -p 10 -f
	echo "${PCODE}/bowtie-0.12.9/bowtie -q -S --best --strata -v 2 -m 1 -p 10 ${PINDEX} ${FILE} ${PDATA}/${BASEFILE}.sam" >>${PDATA}/${SCRIPT};
        
	#echo "${PCODE}/bowtie-0.12.9/bowtie -q -S --best --strata -m 1 -p 15 --chunkmbs 1024 ${PINDEX} ${FILE} ${PDATA}/${BASEFILE}.sam" >>${PDATA}/${SCRIPT};
	echo "${PCODE}/samtools/samtools view -bS ${PDATA}/${BASEFILE}.sam > ${PDATA}/${BASEFILE}.bam" >>${PDATA}/${SCRIPT};
        echo "${PCODE}/samtools/samtools sort ${PDATA}/${BASEFILE}.bam ${PDATA}/${BASEFILE}.sorted" >>${PDATA}/${SCRIPT};
	echo "${PCODE}/samtools/samtools index ${PDATA}/${BASEFILE}.sorted.bam" >>${PDATA}/${SCRIPT};
       
       	nice -5 qsub -e ${PDATA}/bowtie_se_${ID}.err -o ${PDATA}/bowtie_se_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT}; 
	#find . -empty -type f -print0 | xargs -0 echo rm;
done
