#!/bin/bash

#convert hapmap raw to vcf
#add annotation?


#java -Xmx2g -jar GenomeAnalysisTK.jar \
#  -R ref.fasta \
#  -T VariantsToVCF \
#  -o output.vcf \
#  --variant:RawHapMap input.hapmap \
#  --dbsnp dbsnp.vcf

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input raw hapmap .txt files"
        exit
fi

PDATA=$1;
#------------------------
#hard coded paths - code
#------------------------
PCODE_GATK="/net/isi-scratch/giuseppe/tools/GenomeAnalysisTK-2.5-2";
#------------------------
#hard coded paths - data
#------------------------
#human hg19 genome
#PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";
#PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg18/hg18.fa";
PFASTA_GATK="/net/isi-scratch/giuseppe/GATK_RESOURCES/hg18/Homo_sapiens_assembly18.fasta";
PDATA_GATK="/net/isi-scratch/giuseppe/GATK_RESOURCES/hg18";
#PFASTA_GATK="/net/isi-scratch/giuseppe/GATK_RESOURCES/b36/human_b36_both.fasta";
#PDATA_GATK="/net/isi-scratch/giuseppe/GATK_RESOURCES/b36";

#temp
TMP="/tmp"; #GATK gives "too many open files" exception here if you use too many threads
#key - see http://www.broadinstitute.org/gatk/guide/tagged?tag=phone-home
PKEY="/net/isi-scratch/giuseppe/GATK_RESOURCES/giuseppe.gallone_dpag.ox.ac.uk.key";

for FILE in ${PDATA}/*.txt;
        do
        FILE_ID=`basename ${FILE} ".txt"`;
	#FILE_ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        SCRIPT=${FILE_ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

       #-------------------------------
	#GATK VariantsToVCF
        #-------------------------------
        #removed
        # --dbsnp ${PDATA_GATK}/dbsnp_137.hg18.vcf \\
        echo 'echo "START: raw hapmap to vcf conversion";' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "java -Xmx2g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
              -T VariantsToVCF \\
              -R ${PFASTA_GATK} \\
              --variant:RawHapMap  ${FILE}\\
              -et NO_ET \\
              -K ${PKEY} \\
              -o ${PDATA}/${FILE_ID}.vcf" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'echo "DONE: raw hapmap to vcf conversion";' >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};

        #-----------------------
        #GATK Variant Annotator
        #----------------------
        #removed: --useAllAnnotations requires the snpeff annotator\\
        #-A SnpEff \\
        #--dbsnp ${PDATA_GATK}/dbsnp_137.hg19.vcf \\
        #--alwaysAppendDbsnpId \\
        #echo 'echo "START: GATK Variant Annotator";' >>${PDATA}/${SCRIPT};
        #echo '' >>${PDATA}/${SCRIPT};
        #echo "java -Xmx2g -jar ${PCODE_GATK}/GenomeAnalysisTK.jar \\
        #      -T VariantAnnotator \\
        #      -R ${PFASTA_GATK} \\
        #      --variant ${PDATA}/${FILE_ID}.vcf \\
        #      --dbsnp ${PDATA_GATK}/dbsnp_137.hg18.vcf \\
        #      -A FisherStrand -A HaplotypeScore -A InbreedingCoeff -A Coverage \\
        #      -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth \\
        #      -A RMSMappingQuality -A ReadPosRankSumTest -A SpanningDeletions \\
        #      -A AlleleBalance -A AlleleBalanceBySample \\
        #      --alwaysAppendDbsnpId \\
        #      -nt 6 \\
        #      -o ${PDATA}/${FILE_ID}_annotated.vcf \\
        #      -et NO_ET \\
        #      -K ${PKEY}" >>${PDATA}/${SCRIPT};
        #echo '' >>${PDATA}/${SCRIPT};
        #echo 'echo "DONE: GATK Variant Annotator";' >>${PDATA}/${SCRIPT};
        #echo '' >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/${FILE_ID}.err -o ${PDATA}/${FILE_ID}.out -v "BASH_ENV=~/.bashrc" -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
