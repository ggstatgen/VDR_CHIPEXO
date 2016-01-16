#!/bin/bash
#script to generate 2 bowtie indexes for the personalised genomes created by vcf2diploid
#rather stupid at them moment, you need to run it once for each sample

#sample run
#bowtie-build  1_NA19191_paternal.fa,2_NA19191_paternal.fa,3_NA19191_paternal.fa,4_NA19191_paternal.fa,5_NA19191_paternal.fa,6_NA19191_paternal.fa,7_NA19191_paternal.fa,8_NA19191_paternal.fa,9_NA19191_paternal.fa,10_NA19191_paternal.fa,11_NA19191_paternal.fa,12_NA19191_paternal.fa,13_NA19191_paternal.fa,14_NA19191_paternal.fa,15_NA19191_paternal.fa,16_NA19191_paternal.fa,17_NA19191_paternal.fa,18_NA19191_paternal.fa,19_NA19191_paternal.fa,20_NA19191_paternal.fa,21_NA19191_paternal.fa,22_NA19191_paternal.fa PatRef &> bowtie_build_pat.txt &


PCODE="/net/isi-scratch/giuseppe/tools";


if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH_DATA> <SAMPLE_ID>"
        echo "PATH_DATA - directory containing the paternal and maternal fasta files for the sample, generated with do_vcf2diploid_cluster.sh"
	echo "ID of the LCL sample"
        exit
fi

PDATA=$1;
ID=$2;
SEP=",";

#for i in `ls -v *.fas`; do echo $i; done;

for FILE in `ls -v ${PDATA}/*paternal.fa`;
        do
        INPUTPAT+=${FILE}${SEP};
done
INPUTPAT=${INPUTPAT%?};
SCRIPTPAT=bowtie_idx_pat_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPTPAT};
echo '' >>${PDATA}/${SCRIPTPAT};
echo 'source activate' >>${PDATA}/${SCRIPTPAT};
echo '' >>${PDATA}/${SCRIPTPAT};
echo "${PCODE}/bowtie-0.12.9/bowtie-build ${INPUTPAT} ${PDATA}/PatRef" >> ${PDATA}/${SCRIPTPAT};
nice -5 qsub -e ${PDATA}/bowtie_idx_pat_${ID}.err -o ${PDATA}/bowtie_idx_pat_${ID}.out -q newnodes.q ${PDATA}/${SCRIPTPAT};
rm ${PDATA}/${SCRIPTPAT};


for FILE in `ls -v ${PDATA}/*maternal.fa`;
        do
        INPUTMAT+=${FILE}${SEP};
done
INPUTMAT=${INPUTMAT%?};
SCRIPTMAT=bowtie_idx_mat_${ID}.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPTMAT};
echo '' >>${PDATA}/${SCRIPTMAT};
echo 'source activate' >>${PDATA}/${SCRIPTMAT};
echo '' >>${PDATA}/${SCRIPTMAT};
echo "${PCODE}/bowtie-0.12.9/bowtie-build ${INPUTMAT} ${PDATA}/MatRef" >> ${PDATA}/${SCRIPTMAT};
nice -5 qsub -e ${PDATA}/bowtie_idx_mat_${ID}.err -o ${PDATA}/bowtie_idx_mat_${ID}.out -q newnodes.q ${PDATA}/${SCRIPTMAT};
rm ${PDATA}/${SCRIPTMAT};
