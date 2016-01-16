#!/bin/bash
#part of a series of script to do snp calling from Chip-exo data
#input: BAM files produced by the phase 1 GATK pipeline
#output: vcf raw and filtered SNP calls


############
#TASKS
############
#6. GATK - Unified Genotyper
#7. GATK - VCF filter

#command line arguments:
#1 data path (no final backslash)
#2 Indelrealigner Smith-Waterman flag

#The scripts uses the bam files produced by pipeline 1 (sorted bams, indexes)
#the script uses the bams' base filename as the seed name for all its output
#BWA and BWA+STAMPY have been tested. All the .bam should be in the same directory
#all output files are currently stored in the INPUT directory


if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH_DATA>  <METHOD> <INTERVAL_FILE>"
        echo "PATH_DATA - directory containing the input .bam files"
        echo "METHOD   - Whether to use the UnifiedGenotyper (0) or the HaplotypeCaller (1) for Variant Calling"
        echo "INTERVAL_FILE - name of the interval file for GATK -L option"
        exit
fi

PDATA=$1;
METHOD=$2;
INTERVAL=$3;

if [ ${METHOD} == 0 ]; 
then
	MODEL="USE_READS";
else
	MODEL="USE_SW";
fi

#------------------------
#hard coded paths - code
#------------------------
PCODE_GATK="/net/isi-scratch/giuseppe/tools/GenomeAnalysisTK-2.5-2";
#PCODE_ANNO="/net/isi-scratch/giulio/ANNOVAR";

#------------------------
#hard coded paths - data
#------------------------
PDATA_GATK="/net/isi-scratch/giuseppe/GATK_RESOURCES/hg19";
PDATA_ANNO="/net/isi-scratch/giulio/Tools/ANNOVAR";
#human hg19 genome
PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";
#temp
TMP="/tmp"; #must create subdirs here - GATK gives "too many open files" exception
#key - see http://www.broadinstitute.org/gatk/guide/tagged?tag=phone-home
PKEY="/net/isi-scratch/giuseppe/GATK_RESOURCES/giuseppe.gallone_dpag.ox.ac.uk.key";



for FILE in ${PDATA}/*.bam;
        do
        FILENAME=`basename ${FILE} ".bam"`;
	FILE_ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        SCRIPT=GATK_PH2_${FILE_ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};


	if [ ${METHOD} == 0 ];
	then

        #-----------------------
        #GATK Unified genotyper
        #----------------------
        #this function accepts both -nct and -nt
        #removed:
        #--output_mode EMIT_ALL_SITES \\
        #--genotyping_mode DISCOVERY \\
        echo 'echo "START: GATK Unified genotyper";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx16g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T UnifiedGenotyper \\
              -glm BOTH \\
              -R ${PFASTA} \\
              -L ${INTERVAL} \\
              -I ${FILE} \\
              -D ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
              -o ${PDATA}/${FILENAME}.vcf \\
              -nt 6\\
              -nct 6\\
              -et NO_ET \\
              -K ${PKEY} \\
              -stand_call_conf 10 \\
              -stand_emit_conf 10 \\
              -metrics ${PDATA}/${FILENAME}.metrics \\
              -A Coverage \\
              -A AlleleBalance \\
              -A AlleleBalanceBySample" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK Unified genotyper";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	
	else

        #-----------------------
        #GATK Haplotype Caller
        #----------------------
	#removed:
	#--output_mode EMIT_ALL_SITES \\
	#--genotyping_mode DISCOVERY \\
        echo 'echo "START: GATK Haplotype Caller";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx16g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T HaplotypeCaller \\
              -R ${PFASTA} \\
              -L ${INTERVAL} \\
              -I ${FILE} \\
              -D ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
              -o ${PDATA}/${FILENAME}.vcf \\
              -stand_call_conf 10 \\
              -stand_emit_conf 10 \\
              -et NO_ET \\
              -K ${PKEY}" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: GATK Haplotype Caller";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	fi

	#-----------------------
	#GATK VCF Filter
	#-----------------------
	#echo 'echo "START: GATK VCF Filter";' >>${PDATA}/${SCRIPT};
	#echo '' >>${PDATA}/${SCRIPT};
	#echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
        #      -R ${PFASTA} \\
        #      -et NO_ET \\
        #      -K ${PKEY} \\
        #      -T VariantFiltration \\
        #      --variant ${PDATA}/${FILE_ID}.vcf \\
        #      -o ${PDATA}/${FILE_ID}_filtered.vcf \\
        #      --clusterWindowSize 10 \\" >>${PDATA}/${SCRIPT}; 
	#echo '              --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
        #      --filterName "HARD_TO_VALIDATE" \
        #      --filterExpression "MQ <40" \
        #      --filterName "Map_qual" \
        #      --filterExpression "DP < 10" \
        #      --filterName "LowCoverage" \
        #      --filterExpression "QUAL < 30.0" \
        #      --filterName "VeryLowQual" \
        #      --filterExpression "QUAL > 30.0 && QUAL < 50.0" \
        #      --filterName "LowQual" \
        #      --filterExpression "QD < 2.0" \
        #      --filterName "LowQD" \
        #      --filterExpression "FS > 60" \
        #      --filterName "StrandBias_snp" \
        #      --filterExpression "FS > 200" \
        #      --filterName "StrandBias_indel" \
        #      --filterExpression "HaplotypeScore > 13.0" \
        #      --filterName "Haplotype_score" \
        #      --filterExpression "MQRankSum < -12.5" \
        #      --filterName "MQRanksum" \
        #      --filterExpression "ReadPosRankSum < -8.0" \
        #      --filterName "ReadPosRankSum_snp" \
        #      --filterExpression "ReadPosRankSum < -20.0" \
        #      --filterName "ReadPosRankSum_indel" ' >>${PDATA}/${SCRIPT};
	#echo '' >>${PDATA}/${SCRIPT};
	#echo 'echo "DONE: GATK VCF Filter";' >>${PDATA}/${SCRIPT};
	#echo '' >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/${FILE_ID}.err -o ${PDATA}/${FILE_ID}.out -v "BASH_ENV=~/.bashrc" -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
