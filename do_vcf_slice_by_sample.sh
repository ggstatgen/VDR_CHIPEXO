#!/bin/bash

#Usually VCF files from 1000genomes come with all info about ALL samples in on .vcf
#This scrips accepts as input one such file and creates one job per sample, to output one vcf per sample

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <VCF.GZ> <OUT>"
        echo "VCF - absolute path to full vcf.gz file to slice by sample"
        echo "OUT - output data path"
        exit
fi

INPUT_VCF=$1;
OUT=$2;
BASENAME=`basename ${INPUT_VCF} ".vcf.gz"`;
#PCODE="/net/isi-scratch/giuseppe/tools/vcftools_0.1.10/bin";
PCODE="/net/isi-scratch/giuseppe/tools/vcftools_0.1.12b/bin";
PCODE_CGAT="/net/isi-cgat/ifs/apps/bio/vcftools-0.1.8a/bin";

#trios and unrelated
#format: maternal, paternal, child (expected by AlleleSeq pipeline)
CEU_13291="NA07045,NA06986,NA06997";
YRI_Y110="NA19214,NA19213,NA19215";
YRI_Y111="NA19190,NA19189,NA19191";
YRI_Y116="NA19235,NA19236,NA19237";
YRI_Y120="NA19247,NA19248,NA19249";
CEU_UNREL="NA10846,NA10847,NA11829,NA11832,NA12383,NA12489,NA11919,NA11918,NA07029,NA12716,NA12264,NA12872,NA12752,NA06989,NA10831";

#28/1/14
#alleleseq: this is the only cell for which I have full sequencing AND parent data (I have full sequencing for some YRI parents, but don't have parents genotypes)
#maternal paternal individual
CEU_1334="NA12239,NA12146,NA10847";

#The following two groups are for the alleleseq results on the 20 samples: 15 used genotypes from 1kg; 5 used genotypes from the IMPUTE2 pipeline
ALLELESEQ_IMPUTE="NA07029,NA10831,NA11832,NA19191,NA19249";
ALLELESEQ_1KG="NA06986,NA06989,NA10847,NA11829,NA11919,NA12383,NA12489,NA12872,NA19189,NA19190,NA19213,NA19235,NA19236,NA19247,NA19248";

#sets
#(divided by dataset, eg, omni 2.5, 1000gen)
#FIFTEEN_1K_GEN="NA06986,NA06989,NA10847,NA11829,NA11919,NA12383,NA12489,NA12716,NA12872,NA19189,NA19190,NA19213,NA19235,NA19236,NA19247,NA19248";
FIFTEEN_1K_GEN="NA06986,NA06989,NA10847,NA11829,NA11919,NA12383,NA12489,NA12872,NA19189,NA19190,NA19213,NA19235,NA19236,NA19247,NA19248";

OMNI_25_GEN="NA06986,NA06989,NA10847,NA11829,NA11832,NA11918,NA11919,NA12383,NA12489,NA12716,NA12872,NA19189,NA19190,NA19191,NA19213,NA19214,NA19215,NA19235,NA19236,NA19237,NA19247,NA19248,NA19249";
ALL="NA06986,NA06989,NA06997,NA07029,NA07045,NA10831,NA10846,NA10847,NA11829,NA11832,NA11918,NA11919,NA12264,NA12383,NA12489,NA12716,NA12752,NA12872,NA19189,NA19190,NA19191,NA19213,NA19214,NA19215,NA19235,NA19236,NA19237,NA19247,NA19248,NA19249";

#for i in NA11832 NA07029 NA10831 NA12752
#for i in NA19191 NA19249
#for i in NA06986 NA06989 NA10847 NA11829 NA11919 NA12383 NA12489 NA12872 NA19189 NA19190 NA19213 NA19235 NA19236 NA19247 NA19248
#for i in NA06986 NA06989 NA06997 NA07029 NA07045 NA10831 NA10846 NA10847 NA11829 NA11832 NA11918 NA11919 NA12264 NA12383 NA12489 NA12716 NA12752 NA12872 NA19189 NA19190 NA19191 NA19213 NA19214 NA19215 NA19235 NA19236 NA19237 NA19247 NA19248 NA19249
#for i in NA06989 NA10847 NA19189 NA19190
#do 
#	SCRIPT=vcft_${i}.sh;
#        echo '#!/bin/bash' >>${OUT}/${SCRIPT};
#        echo '' >>${OUT}/${SCRIPT};
#        echo "${PCODE}/vcf-subset -e -c ${i} ${INPUT_VCF} | bgzip -c > ${OUT}/${BASENAME}_${i}.vcf.gz" >>${OUT}/${SCRIPT};
#	nice -5 qsub -e ${OUT}/vcft_${i}.err -o ${OUT}/vcft_${i}.out -q medium_jobs.q ${OUT}/${SCRIPT};
#	rm ${OUT}/${SCRIPT};
#done

#for i in CEU_13291 YRI_Y110 YRI_Y111 YRI_Y116 YRI_Y120
#for i in CEU_1334
#for i in YRI_Y110 YRI_Y111 YRI_Y116 YRI_Y120
#do
#	SCRIPT=vcft_trios_${i}.sh;
#        echo '#!/bin/bash' >>${OUT}/${SCRIPT};
#        echo '' >>${OUT}/${SCRIPT};
#        echo "${PCODE}/vcf-subset -e -c ${!i} ${INPUT_VCF} | bgzip -c > ${OUT}/${BASENAME}_${i}.vcf.gz" >>${OUT}/${SCRIPT};
#        nice -5 qsub -e ${OUT}/vcft_${i}.err -o ${OUT}/vcft_${i}.out -v "BASH_ENV=/home/giuseppe/.bashrc,PATH=${PATH}:/net/isi-scratch/giuseppe/tools/vcftools_0.1.10/lib/perl5/site_perl:${PCODE}"  -q newnodes.q ${OUT}/${SCRIPT};
#        #rm ${OUT}/${SCRIPT};
#done

#SCRIPT=vcft_ALLELESEQ_IMPUTE.sh;
#echo '#!/bin/bash' >>${OUT}/${SCRIPT};
#echo '' >>${OUT}/${SCRIPT};
#echo "${PCODE}/vcf-subset -e -c ${ALLELESEQ_IMPUTE} ${INPUT_VCF} | bgzip -c > ${OUT}/${BASENAME}_ALLELESEQ_IMPUTE.vcf.gz" >>${OUT}/${SCRIPT};
#nice -5 qsub -e ${OUT}/vcft_ALLELESEQ_IMPUTE.err -o ${OUT}/vcft_ALLELESEQ_IMPUTE.out -q medium_jobs.q ${OUT}/${SCRIPT};
#rm ${OUT}/${SCRIPT};

SCRIPT=vcft_ALLELESEQ_1KG.sh;
echo '#!/bin/bash' >>${OUT}/${SCRIPT};
echo '' >>${OUT}/${SCRIPT};
echo "${PCODE}/vcf-subset -e -c ${ALLELESEQ_1KG} ${INPUT_VCF} | bgzip -c > ${OUT}/${BASENAME}_ALLELESEQ_1KG.vcf.gz" >>${OUT}/${SCRIPT};
nice -5 qsub -e ${OUT}/vcft_ALLELESEQ_1KG.err -o ${OUT}/vcft_ALLELESEQ_1KG.out -q medium_jobs.q ${OUT}/${SCRIPT};
rm ${OUT}/${SCRIPT};


#SCRIPT=vcft_FIFTEEN_1K_GEN.sh;
#echo '#!/bin/bash' >>${OUT}/${SCRIPT};
#echo '' >>${OUT}/${SCRIPT};
#echo "${PCODE}/vcf-subset -e -c ${FIFTEEN_1K_GEN} ${INPUT_VCF} | bgzip -c > ${OUT}/${BASENAME}_FIFTEEN_1K_GEN.vcf.gz" >>${OUT}/${SCRIPT};
#nice -5 qsub -e ${OUT}/vcft_FIFTEEN_1K_GEN.err -o ${OUT}/vcft_FIFTEEN_1K_GEN.out -q medium_jobs.q ${OUT}/${SCRIPT};
#rm ${OUT}/${SCRIPT};

#SCRIPT=vcft_OMNI_25_GEN.sh;
#echo '#!/bin/bash' >>${OUT}/${SCRIPT};
#echo '' >>${OUT}/${SCRIPT};
#echo "${PCODE}/vcf-subset -e -c ${OMNI_25_GEN} ${INPUT_VCF} | bgzip -c > ${OUT}/${BASENAME}_OMNI_25_GEN.vcf.gz" >>${OUT}/${SCRIPT};
#nice -5 qsub -e ${OUT}/vcft_OMNI_25_GEN.err -o ${OUT}/vcft_OMNI_25_GEN.out -q newnodes.q ${OUT}/${SCRIPT};
#rm ${OUT}/${SCRIPT};

#SCRIPT=vcft_YRI_Y116.sh;
#echo '#!/bin/bash' >>${OUT}/${SCRIPT};
#echo '' >>${OUT}/${SCRIPT};
#echo "${PCODE}/vcf-subset -e -c ${YRI_Y116} ${INPUT_VCF} | bgzip -c > ${OUT}/${BASENAME}_YRI_Y116.vcf.gz" >>${OUT}/${SCRIPT};
#nice -5 qsub -e ${OUT}/vcft_YRI_Y116.err -o ${OUT}/vcft_YRI_Y116.out -q newnodes.q ${OUT}/${SCRIPT};
#rm ${OUT}/${SCRIPT};
