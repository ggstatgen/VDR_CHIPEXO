#!/bin/bash

#use this if you already have bam files
PDATA1='/net/isi-scratch/giuseppe/VDR/INITIAL/04_MY_FASTQ_CUTADAPT_6_15_RC/01_FASTQ';
PDATA2='/net/isi-scratch/giuseppe/VDR/RAW/28_SVR11_hg19/01_FASTQ';
PDATA3='/net/isi-scratch/giuseppe/VDR/RAW/38_SVR11_hg18/01_FASTQ';
PDATA4='/net/isi-scratch/giuseppe/VDR/RAW/34_SVR11_hg19_LOW_QUALITY/01_FASTQ';
POUT='/net/isi-scratch/giuseppe/VDR/POOLED_PLUS_34';

for FILE1 in ${PDATA1}/*.fastq;
        do 
        BASENAME=`basename ${FILE1} ".fastq"`;
	ID1=`echo ${FILE1} | egrep -o "GM[0-9]*"`;
	for FILE2 in ${PDATA2}/*.fastq;
		do
		ID2=`echo ${FILE2} | egrep -o "GM[0-9]*"`;
		if [ ${ID2} = ${ID1} ]; then
			for FILE3 in ${PDATA3}/*.fastq;
				do
				ID3=`echo ${FILE3} | egrep -o "GM[0-9]*"`;
				if [ ${ID3} = ${ID2} ]; then
					for FILE4 in ${PDATA4}/*.fastq;
						do
						ID4=`echo ${FILE4} | egrep -o "GM[0-9]*"`;
						if [ ${ID4} = ${ID3} ]; then
							#concatenate the four and write to output
							SCRIPT=${ID1}_${ID2}_${ID3}_${ID4}.sh;
							echo '#!/bin/bash' >>${POUT}/${SCRIPT};
				        		echo '' >>${POUT}/${SCRIPT};
				        		echo "cat ${FILE1} ${FILE2} ${FILE3} ${FILE4} > ${POUT}/${BASENAME}.fastq" >>${POUT}/${SCRIPT};
							nice -5 qsub -e ${POUT}/${ID1}.err -o ${POUT}/${ID1}.out -q newnodes.q ${POUT}/${SCRIPT};
						fi
					done
				fi
			done
		fi
	done

done
