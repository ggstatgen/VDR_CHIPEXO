#!/bin/bash

#script to run the allele sequ pipeline on the cluster, one instance per sample
#you must have unzipped an instance of the alleleseq pipeline in the same dir
#this script creates, FOR EACH SAMPLE: 
#1) one copy of the PIPELINE DIR instance per sample, and customises the PIPELINE.mk file per sample before launching it
#2) ONE RESULT DIR instance per sample: in here, a link to the fastq.gz will be created, and all output will end here

#The script prepares the following per sample- customised makefile:

#PL:=$(BASE)/Pipeline
#SNPS:=$(BASE)/GM12878/snp.calls
#CNVS:=$(BASE)/GM12878/cnv.calls
#BNDS:=hits.bed
#MAPS:=$(BASE)/GM12878/%s_NA12878.map
#FDR_SIMS:=5
#FDR_CUTOFF:=0.1
#
#sourcefiles := $(wildcard *.fastq.gz)
#countfiles := $(subst .fastq.gz,.cnt,$(sourcefiles))
#
#%.cnt:%.fastq.gz
#        bash -c "python $(PL)/MergeBowtie.py \
#           <($(PL)/filter_input.sh $(PL) $< | bowtie --best --strata -v 2 -m 1 -f AltRefFather/AltRefFather - ) \
#           <($(PL)/filter_input.sh $(PL) $< | bowtie --best --strata -v 2 -m 1 -f AltRefMother/AltRefMother - ) \
#           $(MAPS) | python $(PL)/SnpCounts.py $(SNPS) - $(MAPS) $@"
#
#all: interestingHets.txt
#
#check:
#        @echo $(sourcefiles)
#
#counts.txt: $(countfiles)
#        python $(PL)/CombineSnpCounts.py 5 $(SNPS) $(BNDS) $(CNVS) counts.txt counts.log $(countfiles)
#
## calculate false discovery rates
#FDR.txt: counts.txt
#        python $(PL)/FalsePos.py counts.txt $(FDR_SIMS) $(FDR_CUTOFF) > FDR.txt
#
#interestingHets.txt: counts.txt FDR.txt
#        awk -f $(PL)/finalFilter.awk thresh=$(shell awk 'END {print $$6}' FDR.txt) < counts.txt > interestingHets.txt
#
#clean:
#        @rm -f FDR.txt interestingHets.txt counts.txt
#
#cleanall: clean
#        @rm -f *.cnt
#
#.DELETE_ON_ERROR:

#        echo '%.cnt:%.fastq.gz'  >>${PIPELINEDIR}/${PIPELINESCRIPT}
#        echo '  bash -c "python $(PL)/MergeBowtie.py \ '  >>${PIPELINEDIR}/${PIPELINESCRIPT}
#        echo "           <(\$(PL)/filter_input.sh \$(PL) $< | bowtie --best --strata -v 2 -m 1 -f ${DIPLOID_GENOME_PATH}/PatRef - ) \ "  >>${PIPELINEDIR}/${PIPELINESCRIPT};
#        echo "           <(\$(PL)/filter_input.sh \$(PL) $< | bowtie --best --strata -v 2 -m 1 -f ${DIPLOID_GENOME_MATH}/MatRef - ) \ "  >>${PIPELINEDIR}/${PIPELINESCRIPT};
#        echo '          $(MAPS) | python $(PL)/SnpCounts.py $(SNPS) - $(MAPS) $@" '  >>${PIPELINEDIR}/${PIPELINESCRIPT};



if [ ! $# == 5 ]; then
        echo "Usage: `basename $0` <PATH_BASE> <PATH_SNV> <PATH_CNV> <PEAK_FILE> <FDR>"
        echo "<PATH_BASE> absolute path for the OUTPUT directory (where the Alleleseq dir and NAxxxxx.fastq.gz links are and where you want to spawn the output)"
	echo "<PATH_SNV> path to .snv files in alleleseq format "
	echo "<PATH_CNV> path to .cnv files"
	echo "<PEAK_FILE> consensus bed peak file"
	echo "<FDR> FDR cut-off (eg 0.1)"
	echo "NOTE1: personalised genomes are assumed to be under /net/isi-mirror/1kg_alignments/PERSONALISED_GENOMES";
	echo "NOTE2: .fastq.gz reads are assumed to be under /net/isi-scratch/giuseppe/VDR/POOLED_4/OTHER/01_FASTQ";
	echo "If not, modify script";
	echo "NOTE3: the script will attempt to create ONE DIRECTORY per sample by COPYING THE alleleSeq directory"
	echo "SO MAKE SURE IT IS THERE"
        exit
fi

PDATA=$1;
PSNV=$2;
PCNV=$3;
PEAKS=$4;
FDR=$5;

PFASTQ="/net/isi-scratch/giuseppe/VDR/POOLED_4/OTHER/01_FASTQ";
#PFASTQ="/net/isi-mirror/1kg_alignments/RNA";
P_DIP_GEN="/net/isi-mirror/1kg_alignments/PERSONALISED_GENOMES";
BOWTIEP="/net/isi-scratch/giuseppe/tools/bowtie-0.12.9";

for FILE in ${PDATA}/*.fastq.gz;
        do
	ID=`echo ${FILE} | grep -Po "NA\d{5}"`; 
	
	DIPLOID_GENOME_PATH=${P_DIP_GEN}/${ID};
	DIPLOID_GENOME_MAPS=${P_DIP_GEN}/${ID}/%s_${ID}.map;

	#delete previously existing Alleleseq and result directory and create new ones by copying and renaming an existing alleleSeq directory
	PIPELINEDIR="${PDATA}/${ID}_AlleleSeq_pipeline";
	RESULTDIR="${PDATA}/${ID}_AlleleSeq_results";
	if [ -d "$RESULTDIR" ]; then
        	echo "Deleting ${ID}_AlleleSeq_results dir ..";
        	rm -rf $RESULTDIR;
        fi

	if [ -d "$PIPELINEDIR" ]; then
		echo "Deleting ${ID}_AlleleSeq_pipeline dir ..";
	        rm -rf $PIPELINEDIR;
	fi	
        echo "Creating clean ${ID}_AlleleSeq_results dir ..";
	echo "Creating clean ${ID}_AlleleSeq_pipeline dir ..";

	mkdir $RESULTDIR;
	cp -R ${PDATA}/AlleleSeq_pipeline_v1.1 $PIPELINEDIR;

	#link fastq.gz file inside the result dir
	ln -s $FILE $RESULTDIR;

        PIPELINESCRIPT=PIPELINE_${ID}.mk;
	echo "BASE=${RESULTDIR}" >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo "PL:=${PIPELINEDIR}" >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo "BOWTIE:=${BOWTIEP}" >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo "SNPS:=${PSNV}/vcf2snp_${ID}.snv" >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo "CNVS:=${PCNV}/${ID}.cnv" >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo "BNDS:=${PEAKS}" >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo "MAPS:=${DIPLOID_GENOME_MAPS}" >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo "OUT:=${PDATA}" >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo "SAMPLE:=${ID}" >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo 'FDR_SIMS:=5' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo "FDR_CUTOFF:=${FDR}" >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo 'sourcefiles := $(wildcard *.fastq.gz)' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo 'countfiles := $(subst .fastq.gz,.cnt,$(sourcefiles))' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '%.cnt:%.fastq.gz'  >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo "	bash -c \"python \$(PL)/MergeBowtie.py <(\$(PL)/filter_input.sh \$(PL) $< | \$(BOWTIE)/bowtie --best --strata -v 2 -m 1 -f ${DIPLOID_GENOME_PATH}/PatRef - ) <(\$(PL)/filter_input.sh \$(PL) $< | \$(BOWTIE)/bowtie --best --strata -v 2 -m 1 -f ${DIPLOID_GENOME_PATH}/MatRef - ) \$(MAPS) | python \$(PL)/SnpCounts.py \$(SNPS) - \$(MAPS) \$@\"" >>${PIPELINEDIR}/${PIPELINESCRIPT}
	echo '' >>${PIPELINEDIR}/${PIPELINESCRIPT}; 
	echo 'all: interestingHets.txt' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo 'check:' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '	@echo $(sourcefiles)' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo 'counts.txt: $(countfiles)' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '	python $(PL)/CombineSnpCounts.py 5 $(SNPS) $(BNDS) $(CNVS) counts.txt counts.log $(countfiles)' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '# calculate false discovery rates' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo 'FDR.txt: counts.txt' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '	python $(PL)/FalsePos.py counts.txt $(FDR_SIMS) $(FDR_CUTOFF) > FDR.txt' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo 'interestingHets.txt: counts.txt FDR.txt' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '	awk -f $(PL)/finalFilter.awk thresh=$(shell awk '\''END {print $$6}'\'' FDR.txt) < counts.txt > $(OUT)/interestingHets_$(SAMPLE).txt' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo 'clean:' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '	@rm -f FDR.txt interestingHets.txt counts.txt' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo 'cleanall: clean' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '	@rm -f *.cnt' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '' >>${PIPELINEDIR}/${PIPELINESCRIPT};
	echo '.DELETE_ON_ERROR:' >>${PIPELINEDIR}/${PIPELINESCRIPT};


	#now put this in a SHELL SCRIPT WITH SOURCE ACTIVATE And make  -c -v etc etc and launch
	BASHSCRIPT=PIPELINE_${ID}.sh;
        echo '#!/bin/bash' >>${RESULTDIR}/${BASHSCRIPT};
        echo '' >>${RESULTDIR}/${BASHSCRIPT};
        echo 'source activate' >> ${RESULTDIR}/${BASHSCRIPT};
	echo "make --directory ${RESULTDIR} -f ${PIPELINEDIR}/${PIPELINESCRIPT}"  >>${RESULTDIR}/${BASHSCRIPT};

	nice -5 qsub -e ${PDATA}/asb_${ID}.err -o ${PDATA}/asb_${ID}.out -q newnodes.q ${RESULTDIR}/${BASHSCRIPT};
        #rm ${PIPELINEDIR}/${PIPELINESCRIPT};  
done
