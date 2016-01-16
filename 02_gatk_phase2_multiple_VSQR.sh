#!/bin/bash
#part of a series of script to do snp calling from Chip-exo data
#input: BAM files produced by the phase 1 GATK pipeline
#output: vcf raw and filtered SNP calls

#this only sends ONE job to the queue: one call to the SNP caller using ALL the sample BAMS
#-I sample1.bam [-I sample2.bam ...] \

############
#TASKS
############
#1. GATK - Unified Genotyper
#OR
#1. GATK - Haplotype Caller

#2. snpeff annotator
#3. Variant Annotator
#4. VQSR
#5. Apply Recalibration
#6. Variant Eval
#7. vcf to bed

#vsqr workflow
# snp.model <- BuildErrorModelWithVQSR(raw.vcf, SNP)
# indel.model <- BuildErrorModelWithVQSR(raw.vcf, INDEL)
# recalibratedSNPs.rawIndels.vcf <- ApplyRecalibration(raw.vcf, snp.model, SNP)
# analysisReady.vcf <- ApplyRecalibration(recalibratedSNPs.rawIndels.vcf, indel.model, INDEL)

#The scripts uses the bam files produced by pipeline 1 (sorted bams, indexes)
#the script uses the bams' base filename as the seed name for all its output
#BWA and BWA+STAMPY have been tested. All the .bam should be in the same directory
#all output files are currently stored in the INPUT directory

if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <PATH_DATA>  <METHOD> <INTERVAL_FILE> <QUEUE>"
        echo "PATH_DATA - directory containing the input .bam files"
        echo "METHOD   - Whether to use the UnifiedGenotyper (0) or the HaplotypeCaller (1) for Variant Calling"
        echo "INTERVAL_FILE - name of the interval file for GATK -L option"
	echo "QUEUE - which queue to use [medium_jobs|newnodes|fgu217]"
        exit
fi
PDATA=$1;
METHOD=$2;
INTERVAL=$3;
QUEUE=$4;
#------------------------
#hard coded paths
#------------------------
PCODE_GATK="/net/isi-scratch/giuseppe/tools/GenomeAnalysisTK-2.5-2";
PCODE_SNPEFF="/net/isi-scratch/giuseppe/tools/snpEff_2_0_5";
PCODE_VCFTOOLS="/net/isi-scratch/giuseppe/tools/vcftools_0.1.10/bin";
PCODE_ANNO="/net/isi-scratch/giuseppe/tools/annovar";
PDATA_ANNO="/net/isi-scratch/giulio/Tools/ANNOVAR/humandb";
PDATA_GATK="/net/isi-scratch/giuseppe/GATK_RESOURCES/hg19";
#human hg19 genome
PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";
#temp
TMP="/tmp"; #must create subdirs here - GATK gives "too many open files" exception
#key - see http://www.broadinstitute.org/gatk/guide/tagged?tag=phone-home
PKEY="/net/isi-scratch/giuseppe/GATK_RESOURCES/giuseppe.gallone_dpag.ox.ac.uk.key";


for FILE in ${PDATA}/*.bam;
	do
	OPT="-I ";
	SPACE="  ";
	INPUT+=${OPT}${FILE}${SPACE};
	#need to build list of files
	#-I file1 -I file2 -I file3
done

NAME=`basename ${PDATA}`_VQSR;
SCRIPT="GATK_PH2_${NAME}.sh";
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo '' >>${PDATA}/${SCRIPT};

if [ ${METHOD} == 0 ];
then
	BASENAME=${NAME}_UG;
        #-----------------------
        #GATK Unified genotyper
        #----------------------
        #this function accepts both -nct and -nt
	#removed:
	#--genotyping_mode DISCOVERY \\
	#--output_mode EMIT_ALL_SITES \\
        echo 'echo "START: GATK Unified genotyper";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx16g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T UnifiedGenotyper \\
              -glm BOTH \\
              -R ${PFASTA} \\
              -L ${INTERVAL} \\
              ${INPUT} \\
              -D ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
              -o ${PDATA}/${BASENAME}.vcf \\
              -nt 6\\
              -nct 6\\
              -et NO_ET \\
              -K ${PKEY} \\
              --standard_min_confidence_threshold_for_calling 20 \\
              --standard_min_confidence_threshold_for_emitting 20 \\
              -metrics ${PDATA}/${BASENAME}.vcf.metrics \\
              -A FisherStrand -A HaplotypeScore -A InbreedingCoeff -A Coverage \\
              -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth \\
              -A RMSMappingQuality -A ReadPosRankSumTest -A SpanningDeletions \\
              -A AlleleBalance -A AlleleBalanceBySample" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK Unified genotyper";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

else
	BASENAME=${NAME}_HC;
        #-----------------------
        #GATK Haplotype Caller
        #----------------------
        #add annotations
        echo 'echo "START: GATK Haplotype Caller";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx16g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T HaplotypeCaller \\
              -R ${PFASTA} \\
              -L ${INTERVAL} \\
              --standard_min_confidence_threshold_for_calling 20 \\
              --standard_min_confidence_threshold_for_emitting 20 \\
              ${INPUT} \\
              -D ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
              -o ${PDATA}/${BASENAME}.vcf \\
              -A FisherStrand -A HaplotypeScore -A InbreedingCoeff -A Coverage \\
              -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth \\
              -A RMSMappingQuality -A ReadPosRankSumTest -A SpanningDeletions \\
              -A AlleleBalance -A AlleleBalanceBySample \\
              -et NO_ET \\
              -K ${PKEY}" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK Haplotype Caller";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

fi

        #----------------------------
        #snpeff - variant annotation
        #----------------------------
        #creates a vcf which you need to pass to variantannotator
        echo 'echo "START: snpeff 2.0.5 variant annotation";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx4G -jar ${PCODE_SNPEFF}/snpEff.jar \\
              -v -onlyCoding true \\
              -i vcf -o vcf \\
              -c ${PCODE_SNPEFF}/snpEff.config \\
              -s ${PDATA}/${BASENAME}_snpEff_summary.html \\
              hg19 \\
              ${PDATA}/${BASENAME}.vcf > ${PDATA}/${BASENAME}_snpEff_output.vcf" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: snpeff 2.0.5 variant annotation";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        #-----------------------
        #GATK Variant Annotator
        #----------------------
        #removed: --useAllAnnotations requires the snpeff annotator\\
        #-A SnpEff \\
        #--dbsnp ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
        #--alwaysAppendDbsnpId \\
        #For this to work, the input vcf needs to be annotated with the corresponding values (QD, FS, DP, etc.).
        #If any of these values are somehow missing, then VariantAnnotator needs to be run first so that VariantRecalibration can run properly.
        echo 'echo "START: GATK Variant Annotator";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx16g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T VariantAnnotator \\
              -R ${PFASTA} \\
              --variant ${PDATA}/${BASENAME}.vcf \\
              -L ${INTERVAL} \\
              -A SnpEff \\
              --snpEffFile ${PDATA}/${BASENAME}_snpEff_output.vcf \\
              --dbsnp ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
              --alwaysAppendDbsnpId \\
              -nt 6 \\
              -o ${PDATA}/${BASENAME}_snpeff.vcf \\
              -et NO_ET \\
              -K ${PKEY}" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK Variant Annotator";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        #---------------------------
        #GATK Variant Recalibrator - SNP
        #---------------------------
        #[SPECIFY TRUTH AND TRAINING SETS] 
        #[SPECIFY WHICH ANNOTATIONS TO USE IN MODELING] 
        #[SPECIFY WHICH CLASS OF VARIATION TO MODEL] 
	#removed -percentBad 0.01 -minNumBad 1000 \\
        #removed -qual 10.0 \\
	#removed -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP -an HaplotypeScore -an MQ0 \\
        echo 'echo "START: GATK Variant Recalibrator";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx16g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T VariantRecalibrator \\
              -R ${PFASTA} \\
              -L ${INTERVAL}  \\
              -input ${PDATA}/${BASENAME}.vcf \\
              -recalFile ${PDATA}/${BASENAME}.vcf.recal \\
              -tranchesFile ${PDATA}/${BASENAME}.vcf.tranches \\
              --maxGaussians 2 -percentBad 0.05 \\
              -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${PDATA_GATK}/hapmap_3.3.hg19.vcf \\
              -resource:omni,known=false,training=true,truth=true,prior=12.0 ${PDATA_GATK}/1000G_omni2.5.hg19.vcf  \\
              -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${PDATA_GATK}/1000G_phase1.snps.high_confidence.hg19.vcf \\
              -resource:dbsnp,known=true,training=false,truth=false,prior=2.0  ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
              -rscriptFile ${PDATA}/${BASENAME}_plots.R \\
              -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP -an HaplotypeScore \\
              -mode SNP \\
              -nt 6 \\
              -et NO_ET \\
              -K ${PKEY}" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK VariantRecalibrator";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        #-----------------------
        #GATK Apply Recalibration - SNP
        #----------------------
        #   [SPECIFY THE DESIRED LEVEL OF SENSITIVITY TO TRUTH SITES]
        #   [SPECIFY WHICH CLASS OF VARIATION WAS MODELED]
        echo 'echo "START: GATK Apply Recalibration";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx16g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T ApplyRecalibration \\
              -R ${PFASTA} \\
              -L ${INTERVAL} \\
              -input ${PDATA}/${BASENAME}.vcf \\
              -tranchesFile ${PDATA}/${BASENAME}.vcf.tranches \\
              -recalFile ${PDATA}/${BASENAME}.vcf.recal \\
              -o ${PDATA}/${BASENAME}_recal_filtered.vcf \\
              --ts_filter_level 99.9 \\
              -mode SNP \\
              -et NO_ET \\
              -K ${PKEY}" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK Apply Recalibration";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};


        #-----------------------
        #GATK Variant Eval
        #----------------------
        #see http://gatkforums.broadinstitute.org/discussion/48/using-varianteval
        #removed -noST -noEV -EV CountVariants,TiTvVariantEvaluator -ST CompRod,EvalRod,Sample,Novelty,Filter \\
        echo 'echo "START: GATK Variant Eval";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx16g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T VariantEval \\
              -R ${PFASTA} \\
              -L ${INTERVAL} \\
              -nt 6 \\
              -eval ${PDATA}/${BASENAME}_recal_filtered.vcf \\
              -o ${PDATA}/${BASENAME}.gatkreport \\
              --dbsnp ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
              -et NO_ET \\
              -K ${PKEY}" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK Variant Eval";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};


        #-----------------------
        #GATK Produce BEAGLE input
        #----------------------
        #see http://gatkforums.broadinstitute.org/discussion/48/using-varianteval
        #removed -noST -noEV -EV CountVariants,TiTvVariantEvaluator -ST CompRod,EvalRod,Sample,Novelty,Filter \\
        echo 'echo "START: GATK Produce BEAGLE Input";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx2g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T ProduceBeagleInput \\
              -R ${PFASTA} \\
              -L ${INTERVAL} \\
              -V ${PDATA}/${BASENAME}_recal_filtered.vcf \\
              -o ${PDATA}/${BASENAME}.beagle.output \\
              -et NO_ET \\
              -K ${PKEY}" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK Produce BEAGLE input";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        #-----------------------
        #GATK VCF to BED
        #----------------------
        #echo 'echo "START: GATK VCF to BED";' >>${PDATA}/${SCRIPT};
        #echo '' >>${PDATA}/${SCRIPT};
        #echo "java -Xmx16g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
        #      -T VariantsToBinaryPed \\
        #      -R ${PFASTA} \\
        #      -L ${INTERVAL} \\
        #      --variant ${PDATA}/${BASENAME}_hf.vcf \\
        #      --bed ${PDATA}/${BASENAME}_hf.bed \\
        #      --bim ${PDATA}/${BASENAME}_hf.bim \\
        #      --fam ${PDATA}/${BASENAME}_hf.fam \\
        #      --metaData ${PDATA}/${BASENAME}_hf.metaData \\
        #      --minGenotypeQuality 0\\
        #     -et NO_ET \\
        #      -K ${PKEY}" >>${PDATA}/${SCRIPT};
        #echo '' >>${PDATA}/${SCRIPT};
        #echo 'echo "DONE: GATK VCF to BED";' >>${PDATA}/${SCRIPT};
        #echo '' >>${PDATA}/${SCRIPT};

nice -5 qsub -e ${PDATA}/${BASENAME}.err -o ${PDATA}/${BASENAME}.out -v "BASH_ENV=/home/giuseppe/.bashrc,PATH=${PATH}:/net/isi-cgat/ifs/apps/apps/R-2.14.1/bin/Rscript" -q ${QUEUE}.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};
