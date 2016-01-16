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
#1. GATK - VCF filter

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

NAME=`basename ${PDATA}`_HF;
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

#HARD FILTERING
#---
#snp
#---
#QD < 2.0 x
#MQ < 40.0 x
#FS > 60.0 x
#HaplotypeScore > 13.0  x
#MQRankSum < -12.5 x
#ReadPosRankSum < -8.0 x
#------------------
#indels
#------------------
#QD < 2.0
#ReadPosRankSum < -20.0
#InbreedingCoeff < -0.8
#FS > 200.0
#the InbreedingCoeff statistic is a population-level calculation that is only available with 10 or more samples. If you have fewer samples you will need to omit that particular filter statement.
        #-----------------------
        #GATK VCF Filter
        #-----------------------
        echo 'echo "START: GATK VCF Filter";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -R ${PFASTA} \\
              -et NO_ET \\
              -K ${PKEY} \\
              -T VariantFiltration \\
              --variant ${PDATA}/${BASENAME}.vcf \\
              -o ${PDATA}/${BASENAME}_hf.vcf" >>${PDATA}/${SCRIPT}; 
        echo '              --filterExpression "MQ <40" \
              --filterName "Map_qual" \
              --filterExpression "DP < 10" \
              --filterName "LowCoverage" \
              --filterExpression "QUAL < 30.0" \
              --filterName "VeryLowQual" \
              --filterExpression "QUAL > 30.0 && QUAL < 50.0" \
              --filterName "LowQual" \
              --filterExpression "QD < 2.0" \
              --filterName "LowQD" \
              --filterExpression "FS > 60" \
              --filterName "StrandBias_snp" \
              --filterExpression "HaplotypeScore > 13.0" \
              --filterName "Haplotype_score" \
              --filterExpression "MQRankSum < -12.5" \
              --filterName "MQRanksum" \
              --filterExpression "ReadPosRankSum < -8.0" \
              --filterName "ReadPosRankSum_snp" ' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK VCF Filter";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

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
              ${PDATA}/${BASENAME}_hf.vcf > ${PDATA}/${BASENAME}_snpEff_output.vcf" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: snpeff 2.0.5 variant annotation";' >>${PDATA}/${SCRIPT};
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
              -eval ${PDATA}/${BASENAME}_hf.vcf \\
              -o ${PDATA}/${BASENAME}.hf.gatkreport \\
              --dbsnp ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
              -et NO_ET \\
              -K ${PKEY}" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK Variant Eval";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        #-----------------------
        #GATK Select Variants
        #----------------------
        #see http://gatkforums.broadinstitute.org/discussion/48/using-varianteval
        echo 'echo "START: GATK Select Variants";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx2g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T SelectVariants \\
              -R ${PFASTA} \\
              -nt 6 \\
              -V ${PDATA}/${BASENAME}_hf.vcf \\
              -o ${PDATA}/${BASENAME}.hf_PASS.vcf \\
              -ef -env \\
              -et NO_ET \\
              -K ${PKEY}" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK Select Variants";' >>${PDATA}/${SCRIPT};
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
              -V ${PDATA}/${BASENAME}_hf.vcf \\
              -o ${PDATA}/${BASENAME}.beagle.output \\
              -et NO_ET \\
              -K ${PKEY}" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK Produce BEAGLE input";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};


 
nice -5 qsub -e ${PDATA}/${BASENAME}.err -o ${PDATA}/${BASENAME}.out -v "BASH_ENV=/home/giuseppe/.bashrc,PATH=${PATH}:/net/isi-cgat/ifs/apps/apps/R-2.14.1/bin/Rscript" -q ${QUEUE}.q ${PDATA}/${SCRIPT};
rm ${PDATA}/${SCRIPT};
