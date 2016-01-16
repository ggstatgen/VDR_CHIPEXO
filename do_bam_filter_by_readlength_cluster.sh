#!/bin/bash

#script to filter bams to retain only reads longer than a certain threshold


if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <MIN_RL>"
        echo "<PATH> full path for the BAM files"
	echo "<MIN_RL> minimum read length to retain, eg 20"
        exit 
fi

PDATA=$1;
RL=$2;
SAMTOOLS="/home/giuseppe/local/bin/samtools";

for FILE in ${PDATA}/*.bam;
	do ID=`basename ${FILE} ".bam"`;

        SCRIPT=filter_rl_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	
	#extract header
	echo "${SAMTOOLS} view -H ${FILE} > ${PDATA}/${ID}_minRL_${RL}.sam" >> ${PDATA}/${SCRIPT};
	#filter reads
	#you can change the comparison in here to pick ==, <=. ecc
	echo "${SAMTOOLS} view ${FILE} | perl -ne '\$rl_field = (split /\t/)[9]; print \$_ if(length(\$rl_field) >= ${RL});' >> ${PDATA}/${ID}_minRL_${RL}.sam" >>${PDATA}/${SCRIPT};	
	#create bam,sort,ndex
	echo "${SAMTOOLS} view -bS ${PDATA}/${ID}_minRL_${RL}.sam > ${PDATA}/${ID}_minRL_${RL}.bam" >>${PDATA}/${SCRIPT};
	echo "${SAMTOOLS} sort ${PDATA}/${ID}_minRL_${RL}.bam ${PDATA}/${ID}_minRL_${RL}.sorted" >>${PDATA}/${SCRIPT};
	echo "${SAMTOOLS} index ${PDATA}/${ID}_minRL_${RL}.sorted.bam" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/filter_rl_${ID}.err -o ${PDATA}/filter_rl_${ID}.out -q fgu217.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
