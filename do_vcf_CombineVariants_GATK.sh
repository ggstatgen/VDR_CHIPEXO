#!/bin/bash

#use gatk to combine variants obtained from several sources:
#1000 genomes
#omni 2.5 1000 genomes
#hapmap

#variants should all be v37
#just a wrapper for the actual function currently

#walker format:
# java -Xmx2g -jar GenomeAnalysisTK.jar \
#   -R ref.fasta \
#   -T CombineVariants \
#   --variant:foo input1.vcf \
#   --variant:bar input2.vcf \
#   -o output.vcf \
#   -genotypeMergeOptions PRIORITIZE
#   -priority foo,bar


if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH_DATA> <INTERVAL_FILE> <QUEUE>"
        echo "PATH_DATA - directory containing the input .bam files"
        echo "INTERVAL_FILE - name of the interval file for GATK -L option"
	echo "QUEUE - which queue to use [medium_jobs|newnodes|fgu217]"
        exit
fi
PDATA=$1;
INTERVAL=$2;
QUEUE=$3;
#------------------------
#hard coded paths
#------------------------
PCODE_GATK="/net/isi-scratch/giuseppe/tools/GenomeAnalysisTK-2.5-2";
#PCODE_VCFTOOLS="/net/isi-scratch/giuseppe/tools/vcftools_0.1.10/bin";
#PDATA_GATK="/net/isi-scratch/giuseppe/GATK_RESOURCES/hg19";
#human hg19 genome
#PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";
PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/g1k_v37/human_g1k_v37.fasta";
#temp
TMP="/tmp"; #must create subdirs here - GATK gives "too many open files" exception
#key - see http://www.broadinstitute.org/gatk/guide/tagged?tag=phone-home
PKEY="/net/isi-scratch/giuseppe/GATK_RESOURCES/giuseppe.gallone_dpag.ox.ac.uk.key";

SCRIPT="GATK_combine.sh";
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};

        echo "java -Xmx2g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T CombineVariants \\
              -R ${PFASTA} \\
              -L ${INTERVAL} \\
              -nt 6 \\
              --variant:g1k $PDATA/1000g/ALL.phase1_release_v3.20101123.snps_indels_svs.genotypes_FIFTEEN_1K_GEN.vcf.gz \\
              --variant:omni25 $PDATA/omni_25_b37/Omni25_genotypes_2141_samples.b37_OMNI_25_GEN_PASS.vcf.gz \\
              --variant:hapmap27 $PDATA/hapmap_liftover_b37/genotypes_CEUYRI_30_r27_nr.b37_fwd.vcf.gz \\
              -o $PDATA/COMBO.vcf \\
              -genotypeMergeOptions PRIORITIZE \\
              -priority g1k,omni25,hapmap27 \\
              -et NO_ET \\
              -K ${PKEY}" >>${PDATA}/${SCRIPT};

nice -5 qsub -e ${PDATA}/${SCRIPT}.err -o ${PDATA}/${SCRIPT}.out -q ${QUEUE}.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};
