#!/bin/bash


#script to extract bam alignments from sra files
#sra toolkit from here
#http://eutils.ncbi.nih.gov/Traces/sra/?view=software
#docs are here
#http://eutils.ncbi.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std

#typical command
#/net/isi-scratch/giuseppe/tools/sratoolkit.2.3.4-2-centos_linux64/bin/sam-dump --unaligned SRR833574 | samtools view -Sb -o SRR833574.bam -

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "<PATH> data path for the .sra files"
        exit
fi

PDATA=$1;
#PCODE="/net/isi-scratch/giuseppe/tools/sratoolkit.2.3.5-2-centos_linux64/bin";
PCODE="/net/isi-scratch/giuseppe/tools/SRA/linux/gcc/stat/x86_64/rel/bin"
SAMTOOLS="/net/isi-scratch/giuseppe/tools/samtools-0.1.19/samtools";

for FILE in ${PDATA}/*.sra;
	do ID=`basename ${FILE} ".sra"`;

        SCRIPT=sra2bam_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	#echo 'source activate' >>${PDATA}/${SCRIPT};
        echo "${PCODE}/sam-dump --unaligned ${FILE} | ${SAMTOOLS} view -Sb -o ${PDATA}/${ID}.bam -" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/sra2bam_${ID}.err -o ${PDATA}/sra2bam_${ID}.out  -q fgu205.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
