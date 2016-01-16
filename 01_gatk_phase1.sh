#!/bin/bash
#part of a series of script to do snp calling from Chip-exo data
#input: SAM files produced with BWA
#processing via picard includes: bam conversion sorting, duplicate marking and read-group-processing via Picard
#processing via gatk include local realignment around indels, base quality recalibration and removal of reads outside specified interval (output bams will be smaller)
#output: realigned bam,metrics and input for 02_gatk... module, where the actual SNP calling happens

# SAM->BAM->......>BAM ready for SNP calling

#roughly corresponds to the Phase 1 - NGS data processing block under
#http://www.broadinstitute.org/gatk/guide/topic?name=best-practices
#Giuseppe Gallone 3/6/2013

############
#TASKS
############
#1. PICARD - sam -> bam conversion
#2. PICARD - mark PCR duplicates
#3. PICARD - fix Read Group information/generate index

#4. GATK - create targets to realign
#5. GATK - create realigned BAM
#6. GATK - base recalibrator 
#7. GATK - print reads

#command line arguments:
#1 data path (no final backslash)
#2 Indelrealigner Smith-Waterman flag
#3 bed interval file whose reads should be kept. All reads outside these will be removed from the final bams. The interval was obtained with MACS2 using loose thresholds and the union over all cells was obtained

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
PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.84/picard-tools-1.84";
#------------------------
#hard coded paths - data
#------------------------
PDATA_GATK="/net/isi-scratch/giuseppe/GATK_RESOURCES/hg19";
#human hg19 genome
PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";
#temp
TMP="/tmp"; #GATK gives "too many open files" exception here if you use too many threads
#key - see http://www.broadinstitute.org/gatk/guide/tagged?tag=phone-home
PKEY="/net/isi-scratch/giuseppe/GATK_RESOURCES/giuseppe.gallone_dpag.ox.ac.uk.key";

#------------------
#READ GROUP STRINGS
#-------------------
# !!! CHANGE HERE depending on dataset
#
RGID="VDR_sc-1008_LCL" #append to this dynamically from input file 
RGPL="ILLUMINA";  #gatk wants only one of ILLUMINA,SLX,SOLEXA,SOLID,454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN
RGLB="XO111";
RGPU="BWA";  #really don't know what to put in this one
#RGSM=""; #Giulio set this to "Andrew" and "Jean", I cannot get these from the file names, will use some variant of RGID taken from the file name
#
# !!! CHANGE HERE IF NEEDED
#

for FILE in ${PDATA}/*.sam;
        do
        FILENAME=`basename ${FILE} ".sam"`;
	FILE_ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        SCRIPT=GATK_PH1_${FILE_ID}.sh;
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
              OUTPUT=${PDATA}/${FILENAME}.bam \\
              VALIDATION_STRINGENCY=LENIENT" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "*DONE*: Picard SAM to BAM";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        #-----------------------
        #Picard Mark Duplicates
        #-----------------------
        #removed: REMOVE_DUPLICATES=true\\
        echo 'echo "*START*: Picard Mark Duplicates";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx4g -Djava.io.tmpdir=${TMP}\\
              -jar ${PCODE_PICARD}/MarkDuplicates.jar\\
              INPUT=${PDATA}/${FILENAME}.bam\\
              OUTPUT=${PDATA}/${FILENAME}_m.bam\\
              METRICS_FILE=${PDATA}/${FILE_ID}_picard.metrics\\
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
              INPUT=${PDATA}/${FILENAME}_m.bam\\
              OUTPUT=${PDATA}/${FILENAME}_mrg.bam\\
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
	#This determines (small) suspicious intervals which are likely in need of realignment
	#This tool ignores MQ0 reads and reads with consecutive indel operators in the CIGAR string - REMOVE MAPQ=0 reads before doing this?
        #Testing -nt option with n thread and n*2G of memory (Xmx)
	#see http://gatkforums.broadinstitute.org/discussion/1975/recommendations-for-parallelizing-gatk-tools
	echo 'echo "*START*: GATK Create Targets to Realign";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx12g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T RealignerTargetCreator \\
              -R ${PFASTA} \\
              -I ${PDATA}/${FILENAME}_mrg.bam \\
              -L ${INTERVAL}\\
              -nt 8 \\
              -et NO_ET \\
              -K ${PKEY} \\
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
	echo 'echo "*START*: GATK Create Realigned BAM";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
        echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T IndelRealigner \\
              -I ${PDATA}/${FILENAME}_mrg.bam \\
              -R ${PFASTA} \\
              -L ${INTERVAL} \\
              -et NO_ET \\
              -K ${PKEY} \\
              -targetIntervals ${PDATA}/${FILE_ID}_input.bam.list \\
              -model ${MODEL} \\
              -known ${PDATA_GATK}/1000G_phase1.indels.hg19.vcf \\
              -known ${PDATA_GATK}/Mills_and_1000G_gold_standard.indels.hg19.vcf \\
              -o ${PDATA}/${FILE_ID}_mrg_realigned.bam" >>${PDATA}/${SCRIPT}; 
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
              -I ${PDATA}/${FILE_ID}_mrg_realigned.bam\\
              -L ${INTERVAL} \\
              -nct 8\\
              -et NO_ET \\
              -K ${PKEY} \\
              -R ${PFASTA} \\
              -knownSites ${PDATA_GATK}/Mills_and_1000G_gold_standard.indels.hg19.vcf \\
              -knownSites ${PDATA_GATK}/1000G_phase1.indels.hg19.vcf \\
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
              -R ${PFASTA} \\
              -L ${INTERVAL} \\
              -nct 8\\
              -et NO_ET \\
              -K ${PKEY} \\
              -I ${PDATA}/${FILE_ID}_mrg_realigned.bam \\
              -BQSR ${PDATA}/${FILE_ID}_recal_data.grp \\
              -o ${PDATA}/${FILE_ID}_GATK.bam" >>${PDATA}/${SCRIPT}; #marked realigned recalibrated
	echo '' >>${PDATA}/${SCRIPT};
	echo 'echo "*DONE*: GATK Print Reads";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

	#-----------------------
	#GATK Depth of Coverage
	#-----------------------
	#echo 'echo "*START*: GATK Depth of Coverage";' >>${PDATA}/${SCRIPT};
	#echo '' >>${PDATA}/${SCRIPT};
	#echo "java -Xmx4g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
        #      -T DepthOfCoverage \\
        #      -et NO_ET \\
        #      -K ${PKEY} \\
        #      -R ${PFASTA} \\
        #      -L ${INTERVAL} \\
        #      -I ${PDATA}/${FILE_ID}_marked.realigned.recalibrated.bam \\
        #      -o ${PDATA}/${FILE_ID}_coverage" >>${PDATA}/${SCRIPT};
	#echo '' >>${PDATA}/${SCRIPT};
	#echo 'echo "*DONE*: GATK Depth of Coverage";' >>${PDATA}/${SCRIPT};
	#echo '' >>${PDATA}/${SCRIPT};
	
	nice -5 qsub -e ${PDATA}/${FILE_ID}.err -o ${PDATA}/${FILE_ID}.out -v "BASH_ENV=~/.bashrc" -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
