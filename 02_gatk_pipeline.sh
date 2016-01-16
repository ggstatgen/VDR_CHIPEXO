#!/bin/bash
#script to automate GATK analysis of exome sequencing data. 
#current pipeline inspired from
#http://seqanswers.com/wiki/How-to/exome_analysis
#Prerequisite for this pipeline is the picard pipeline in the same directory
#Giuseppe Gallone 5/2/2013


############
#TASKS
############
#1. GATK - create targets to realign
#2. GATK - create realigned BAM
#3. GATK - base recalibrator 
#4. GATK - print reads
#5. GATK - Depth of coverage
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
        echo "Usage: `basename $0` <PATH_DATA>  <SW_FLAG> <INTERVAL_FILE>"
        echo "PATH_DATA - directory containing the input .bam files"
        echo "SW_FLAG   - Indelrealigner flag. Options: 0 (default)/1 (smith-waterman)"
	echo "INTERVAL_FILE - name of the interval file for GATK -L option (required)"
        exit
fi

PDATA=$1;
FLAG=$2;
INTERVAL=$3;

if [ ${FLAG} == 0 ]; 
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
TMP="/tmp";


for FILE in ${PDATA}/*.bam;
        do
        #FILE_ID=`basename ${FILE} ".bam"`;
	FILE_ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        SCRIPT=sample_${FILE_ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        #-------------------------------
	#GATK Create Targets to Realign
        #-------------------------------
	#This determines (small) suspicious intervals which are likely in need of realignment
	#This tool ignores MQ0 reads and reads with consecutive indel operators in the CIGAR string - REMOVE MAPQ=0 reads before doing this!
        #Testing -nt option with n thread and n*2G of memory (Xmx)
	#removed: -L ${PDATA_GATK}/${INTERVAL}\\
	#see http://gatkforums.broadinstitute.org/discussion/1975/recommendations-for-parallelizing-gatk-tools
	echo 'echo "*START*: GATK Create Targets to Realign";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx12g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T RealignerTargetCreator \\
              -R ${PFASTA} \\
              -I ${FILE} \\
              -nt 6\\
              -known ${PDATA_GATK}/1000G_phase1.indels.hg19.vcf \\
              -known ${PDATA_GATK}/Mills_and_1000G_gold_standard.indels.hg19.vcf \\
              -o ${PDATA}/${FILE_ID}_input.bam.list" >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: GATK Create Targets to Realign";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

        #-------------------------------
	#GATK Create Realigned BAM
        #-------------------------------
	#-nt and -nct not possible with this function
	#runs the realigner over the intervals created by RealignerTargetCreator
	#removed: -L ${PDATA_GATK}/${INTERVAL}\\
	echo 'echo "*START*: GATK Create Realigned BAM";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T IndelRealigner \\
              -I ${FILE} \\
              -R ${PFASTA} \\
              -targetIntervals ${PDATA}/${FILE_ID}_input.bam.list \\
              -model ${MODEL} \\
              -known ${PDATA_GATK}/1000G_phase1.indels.hg19.vcf \\
              -known ${PDATA_GATK}/Mills_and_1000G_gold_standard.indels.hg19.vcf \\
              -o ${PDATA}/${FILE_ID}_marked.realigned.bam" >>${PDATA}/${SCRIPT}; 
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: GATK Create Realigned BAM";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	
        #-----------------------
	#GATK Base Recalibrator
	#-----------------------
	#Testing -nct option with 8 threads
        #see http://gatkforums.broadinstitute.org/discussion/1975/recommendations-for-parallelizing-gatk-tools
	#removed: -L ${PDATA_GATK}/${INTERVAL}\\
	echo 'echo "*START*: GATK Base Recalibrator";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T BaseRecalibrator \\
              -I ${PDATA}/${FILE_ID}_marked.realigned.bam\\
              -nct 8\\
              -R ${PFASTA} \\
              -knownSites ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
              -o ${PDATA}/${FILE_ID}_recal_data.grp" >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: GATK Base Recalibrator";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

	#-----------------
	#GATK Print Reads
	#-----------------
	#Testing -nct option with 8 threads
	#removed: -L ${PDATA_GATK}/${INTERVAL}\\
	echo 'echo "*START*: GATK Print Reads";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T PrintReads \\
              -R ${PFASTA} \\
              -nct 8\\
              -I ${PDATA}/${FILE_ID}_marked.realigned.bam \\
              -BQSR ${PDATA}/${FILE_ID}_recal_data.grp \\
              -o ${PDATA}/${FILE_ID}_marked.realigned.recalibrated.bam" >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: GATK Print Reads";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

	#-----------------------
	#GATK Depth of Coverage
	#-----------------------
	#removed: -L ${PDATA_GATK}/${INTERVAL}\\
	echo 'echo "*START*: GATK Depth of Coverage";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T DepthOfCoverage \\
              -R ${PFASTA} \\
              -I ${PDATA}/${FILE_ID}_marked.realigned.recalibrated.bam \\
              -o ${PDATA}/${FILE_ID}_coverage" >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: GATK Depth of Coverage";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

	#-----------------------
	#GATK Unified genotyper
	#----------------------
	#this function accepts both -nct and -nt
	#removed: -L ${PDATA_GATK}/${INTERVAL}\\
	echo 'echo "START: GATK Unified genotyper";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx16g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -glm BOTH \\
              -R ${PFASTA} \\
              -T UnifiedGenotyper \\
              -I ${PDATA}/${FILE_ID}_marked.realigned.recalibrated.bam \\
              -D ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
              -o ${PDATA}/${FILE_ID}.vcf \\
              -nt 6\\
              -nct 6\\
              -metrics ${PDATA}/${FILE_ID}.metrics \\
              -A Coverage \\
              -A AlleleBalance \\
              -A AlleleBalanceBySample" >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "DONE: GATK Unified genotyper";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

	#-----------------------
	#GATK VCF Filter
	#-----------------------
	echo 'echo "START: GATK VCF Filter";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -R ${PFASTA} \\
              -T VariantFiltration \\
              --variant ${PDATA}/${FILE_ID}.vcf \\
              -o ${PDATA}/${FILE_ID}_filtered.vcf \\
              --clusterWindowSize 10 \\" >>${PDATA}/${SCRIPT}; 
	echo '              --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
              --filterName "HARD_TO_VALIDATE" \
              --filterExpression "MQ <40" \
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
              --filterExpression "FS > 200" \
              --filterName "StrandBias_indel" \
              --filterExpression "HaplotypeScore > 13.0" \
              --filterName "Haplotype_score" \
              --filterExpression "MQRankSum < -12.5" \
              --filterName "MQRanksum" \
              --filterExpression "ReadPosRankSum < -8.0" \
              --filterName "ReadPosRankSum_snp" \
              --filterExpression "ReadPosRankSum < -20.0" \
              --filterName "ReadPosRankSum_indel" ' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "DONE: GATK VCF Filter";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

	#-----------------------
	#ANNOVAR 
	#-----------------------
	nice -5 qsub -e ${PDATA}/${FILE_ID}.stderr -o ${PDATA}/${FILE_ID}.stdout -v "BASH_ENV=~/.bashrc" -q newnodes.q ${PDATA}/${SCRIPT};
        #rm ${PDATA}/${SCRIPT};
done

