#!/bin/bash

#use this if you already have bam files
PDATA1='/net/isi-scratch/giuseppe/VDR/04_MY_FASTQ_CUTADAPT_6_15_RC/02_BAM_BWA';
PDATA2='/net/isi-scratch/giuseppe/VDR/RAW/28_SVR11_hg19/03_BAM_BWA_RG';
PDATA3='/net/isi-scratch/giuseppe/VDR/RAW/38_SVR11_hg18/03_BAM_BWA_RG';
POUT='/net/isi-scratch/giuseppe/VDR/RAW/POOLED/02_BAM_BWA';


TMP="/tmp";
PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.84/picard-tools-1.84";

for FILE1 in ${PDATA1}/*.sam;
        do 
        BASENAME=`basename ${FILE1} ".sam"`;
	ID1=`echo ${FILE1} | egrep -o "GM[0-9]*"`;
	for FILE2 in ${PDATA2}/*.sam;
		do
		ID2=`echo ${FILE2} | egrep -o "GM[0-9]*"`;
		if [ ${ID2} = ${ID1} ]; then
			for FILE3 in ${PDATA3}/*.sam;
				do
				ID3=`echo ${FILE3} | egrep -o "GM[0-9]*"`;
				if [ ${ID3} = ${ID2} ]; then
					#run picard merge
					SCRIPT=${ID1}_${ID2}_${ID3}.sh;
					echo '#!/bin/bash' >>${POUT}/${SCRIPT};
				        echo '' >>${POUT}/${SCRIPT};
					echo "java -Xmx4g -Djava.io.tmpdir=${TMP}\\
					-jar ${PCODE_PICARD}/MergeSamFiles.jar\\
					INPUT=${FILE1}\\
					INPUT=${FILE2}\\
					INPUT=${FILE3}\\
					OUTPUT=${POUT}/${BASENAME}_merged.sam\\
					USE_THREADING=true" >>${POUT}/${SCRIPT};
					nice -5 qsub -e ${POUT}/${ID1}.err -o ${POUT}/${ID1}.out -q medium_jobs.q ${POUT}/${SCRIPT};
				fi
			done
		fi
	done

done
