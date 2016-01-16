#!/bin/bash
#script to automate GATK analysis of exome sequencing data. 
#current pipeline inspired from
#http://seqanswers.com/wiki/How-to/exome_analysis
#Giuseppe Gallone 5/2/2013

############
#TASKS
############
#1. PICARD - sam -> bam conversion
#2. PICARD - mark PCR duplicates
#3. PICARD - fix Read Group information
#4. GATK - create targets to realign
#5. GATK - create realigned BAM
#6. GATK - base recalibrator 
#7. GATK - print reads
#8. GATK - Depth of coverage
#9. GATK - Unified Genotyper
#10. GATK - VCF filter

#command line arguments:
#1 data path (no final backslash)
#2 Indelrealigner Smith-Waterman flag

#The script uses the sam files' base filename as the seed name for all of the pipeline outputs
#INPUT: collection of SAM files obtained from an aligner of choice. BWA and BWA+STAMPY have been tested. All the .sam should be in the same directory
#all output files are currently stored in the INPUT directory


if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH_DATA>  <SW_FLAG>"
        echo "PATH_DATA - directory containing the input .sam files"
        echo "SW_FLAG   - Indelrealigner flag. Options: 0 (default)/1 (smith-waterman)"
        exit
fi

PDATA=$1;
FLAG=$2;

if [ ${FLAG} == 0 ]; 
then
	MODEL="USE_READS";
else
	MODEL="USE_SW";
fi

#------------------------
#hard coded paths - code
#------------------------
#PCODE_PICARD="/home/giulio/picard-tools-1.84/picard-tools-1.84";
PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.84/picard-tools-1.84";
PCODE_GATK="/net/isi-scratch/giulio/GATK_first";
PCODE_ANNO="/net/isi-scratch/giulio/ANNOVAR";

#------------------------
#hard coded paths - data
#------------------------
PDATA_GATK="/net/isi-scratch/giulio/GATK_second";
PDATA_ANNO="/net/isi-scratch/giulio/ANNOVAR/humandb";
#human hg19 genome
PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19";
#temp
TMP="/tmp";


#READ GROUP STRINGS
#
# !!! CHANGE HERE IF NEEDED
#
RGID="WTCHG_46599" #append to this dynamically from input file 
RGPL="ILLUMINA";
RGLB="812/12_MPX";
RGPU="STAMPY";  #really don't know what to put in this one
#RGSM=""; #Giulio set this to "Andrew" and "Jean", I cannot get these from the file names, will use some variant of RGID taken from the file name
#
# !!! CHANGE HERE IF NEEDED
#


for FILE in ${PDATA}/*.sam;
        do
        FILE_ID=`basename ${FILE} ".sam"`;

        SCRIPT=sample_${FILE_ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        #----------------
	#Picard SAM->BAM
        #----------------
	echo 'echo "*START*: Picard SAM to BAM";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx4g -Djava.io.tmpdir=${TMP} \\
              -jar ${PCODE_PICARD}/SortSam.jar \\
              SORT_ORDER=coordinate \\
              INPUT=${FILE} \\
              OUTPUT=${PDATA}/${FILE_ID}.bam \\
              VALIDATION_STRINGENCY=LENIENT" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: Picard SAM to BAM";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

        #-----------------------
        #Picard Mark Duplicates
        #-----------------------
	echo 'echo "*START*: Picard Mark Duplicates";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx4g -Djava.io.tmpdir=${TMP}\\
              -jar ${PCODE_PICARD}/MarkDuplicates.jar\\
              INPUT=${PDATA}/${FILE_ID}.bam\\
              OUTPUT=${PDATA}/${FILE_ID}_marked.bam\\
              REMOVE_DUPLICATES=true\\
              METRICS_FILE=${PDATA_GATK}/metrics\\
              VALIDATION_STRINGENCY=LENIENT" >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: Picard Mark Duplicates";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	
	#----------------------
	#Picard Fix Read Group
	#----------------------
        echo 'echo "*START*: Picard Fix Read Group";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx4g -Djava.io.tmpdir=${TMP}\\
              -jar ${PCODE_PICARD}/AddOrReplaceReadGroups.jar\\
              INPUT=${PDATA}/${FILE_ID}_marked.bam\\
              OUTPUT=${PDATA}/${FILE_ID}_marked_rg.bam\\
              RGID=${RGID}_${FILE_ID}\\
              RGLB=${RGLB}\\
              RGPL=${RGPL}\\
              RGPU=${RGPU}\\
              RGSM=sample_${FILE_ID}\\
              CREATE_INDEX=true" >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: Picard Fix Read Group";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	
        #-------------------------------
	#GATK Create Targets to Realign
        #-------------------------------
        #Testing -nt option with n thread and n*2G of memory (Xmx)
	#see http://gatkforums.broadinstitute.org/discussion/1975/recommendations-for-parallelizing-gatk-tools
	echo 'echo "*START*: GATK Create Targets to Realign";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx12g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T RealignerTargetCreator \\
              -R ${PFASTA}/hg19.fa\\
              -L ${PDATA_GATK}/TruSeq_150bp.bed\\
              -I ${PDATA}/${FILE_ID}_marked_rg.bam \\
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
	echo 'echo "*START*: GATK Create Realigned BAM";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T IndelRealigner \\
              -I ${PDATA}/${FILE_ID}_marked_rg.bam \\
              -R ${PFASTA}/hg19.fa \\
              -L ${PDATA_GATK}/TruSeq_150bp.bed\\
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
	echo 'echo "*START*: GATK Base Recalibrator";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T BaseRecalibrator \\
              -I ${PDATA}/${FILE_ID}_marked.realigned.bam\\
              -L ${PDATA_GATK}/TruSeq_150bp.bed\\
              -nct 8\\
              -R ${PFASTA}/hg19.fa\\
              -knownSites ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
              -o ${PDATA}/${FILE_ID}_recal_data.grp" >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: GATK Base Recalibrator";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

	#-----------------
	#GATK Print Reads
	#-----------------
	#Testing -nct option with 8 threads
	echo 'echo "*START*: GATK Print Reads";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T PrintReads \\
              -R ${PFASTA}/hg19.fa \\
              -L ${PDATA_GATK}/TruSeq_150bp.bed\\
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
	echo 'echo "*START*: GATK Depth of Coverage";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T DepthOfCoverage \\
              -R ${PFASTA}/hg19.fa \\
              -L ${PDATA_GATK}/TruSeq_150bp.bed\\
              -I ${PDATA}/${FILE_ID}_marked.realigned.recalibrated.bam \\
              -o ${PDATA}/${FILE_ID}_coverage" >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: GATK Depth of Coverage";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

	#-----------------------
	#GATK Unified genotyper
	#----------------------
	#this function accepts both -nct and -nt
	echo 'echo "START: GATK Unified genotyper";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx16g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -glm BOTH \\
              -R ${PFASTA}/hg19.fa \\
              -T UnifiedGenotyper \\
              -I ${PDATA}/${FILE_ID}_marked.realigned.recalibrated.bam \\
              -D ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
              -o ${PDATA}/${FILE_ID}.vcf \\
              -nt 6\\
              -nct 6\\
              -metrics ${PDATA}/${FILE_ID}.metrics \\
              -L ${PDATA_GATK}/TruSeq_150bp.bed\\
              -A DepthOfCoverage \\
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
              -R ${PFASTA}/hg19.fa \\
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
        rm ${PDATA}/${SCRIPT};
done

