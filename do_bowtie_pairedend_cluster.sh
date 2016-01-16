#!/bin/bash
#script to run bowtie 2.19 with NOTCH data as suggested by Steve Sansom 8/2/2013

#command:
#bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]
#I will use bowtie -1 1.fq -2 2.fq -q --best --strata --m 1 -p 15 -3 2 -5 2 --solexa1.3-quals --chunkmbs 1024
#where
#  -q                 query input files are FASTQ .fq/.fastq (default)
#  --best             hits guaranteed best stratum; ties broken by quality
#  --strata           hits in sub-optimal strata aren't reported (requires --best)
#  -m <int>           suppress all alignments if > <int> exist (def: no limit)
#  -p/--threads <int> number of alignment threads to launch (default: 1)
#  --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 64)
#  -S                 sam output
#  <m1>    Comma-separated list of files containing upstream mates (or the
#          sequences themselves, if -c is set) paired with mates in <m2>
#  <m2>    Comma-separated list of files containing downstream mates (or the
#          sequences themselves if -c is set) paired with mates in <m1>
#  --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3
#  -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads
#  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads


PCODE="/net/isi-scratch/giuseppe/tools";
PINDEX="/net/isi-scratch/giuseppe/indexes/bwa/mm9_no_chrM";
#PINDEX="/net/isi-scratch/giuseppe/indexes/bwa/mm9.fa";
#qui ci vuole il mmus 9, crea un argomento linea di comando


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the fastq.gz files"
        exit
fi

PDATA=$1;

for FILE in ${PDATA}/*.fastq;
	do 
        PAIR=`echo ${FILE} | egrep -o "_2_"`;
        if [ "${PAIR}" = '_2_' ];
        then
                continue;
        fi

        BASEFILE=`basename ${FILE} "_1_sequence.fastq"`; 
	SCRIPT=${BASEFILE}.sh;
 
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#change required options here:
        echo "${PCODE}/bowtie-0.12.9/bowtie -q -S -1 ${PDATA}/${BASEFILE}_1_sequence.fastq -2 ${PDATA}/${BASEFILE}_2_sequence.fastq --best --strata -m 1 -p 15 --chunkmbs 1024 --solexa1.3-quals -3 2 -5 2 ${PINDEX} ${PDATA}/${BASEFILE}.sam" >>${PDATA}/${SCRIPT};
	echo "${PCODE}/samtools/samtools view -bS ${PDATA}/${BASEFILE}.sam > ${PDATA}/${BASEFILE}.bam" >>${PDATA}/${SCRIPT};
        echo "${PCODE}/samtools/samtools sort ${PDATA}/${BASEFILE}.bam ${PDATA}/${BASEFILE}.sorted" >>${PDATA}/${SCRIPT};
	echo "${PCODE}/samtools/samtools index ${PDATA}/${BASEFILE}.sorted.bam" >>${PDATA}/${SCRIPT};
       
       	nice -5 qsub -e ${PDATA}/${BASEFILE}.err -o ${PDATA}/${BASEFILE}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT}; 
	#find . -empty -type f -print0 | xargs -0 echo rm;

done
