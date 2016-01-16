#!/bin/bash

#script to filter SNPS with a given allele count

#Inspired on the rules here
#http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5/READ_ME_beagle_ref
# 1) Markers with <5 copies of the reference allele or <5 copies of the non-reference alleles have been excluded.
#
# 2) Structural variants have been excluded.

#I have excluded structural variants separately

#Here I remove markers with small minor allele frequency
#--maf <float>
#--max-maf <float>
#for example  0.01 < MAF < 0.99

#Include only sites with a Minor Allele Frequency greater than or equal to the "--maf" value and less than or equal to the "--max-maf" value. One of these options may be used without the other. Allele frequency is defined as the number of times an allele appears over all individuals at that site, divided by the total number of non-missing alleles at that site.

#also 
#--hwe <float>


#I want to clean up the 1kg vcfs before getting LD blocks
#method from here
#http://vcftools.sourceforge.net/man_latest.html
#Output a new vcf file from the input vcf file that removes any indel sites
#vcftools --vcf input_file.vcf --remove-indels --recode --recode-INFO-all --out SNPs_only

#To write out the variants that pass through filters use the --recode option. In addition, use --recode-INFO-all to include all data from the INFO fields in the output. By default INFO fields are not written because many filters will alter the variants in a file, rendering the INFO values incorrect.


if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <PATH> <MAF> <MAXMAF> <HWE>"
        echo "PATH - path for the vcf files"
	echo "MAF - minimum minor allele frequency value (eg 0.01) such that  0.01 < MAF < [1-0.01] will be kept";
	echo "MAXMAF - 1-<MAF>";
	echo "HWE - hardy weinberg equilibrium threshold to be used (0.001 used for D' blocks)";
        exit
fi

PDATA=$1;
MAF=$2;
MAXMAF=$3;
HWE=$4;

#location of vcf-concat executable, change if needed
PCODE="/home/giuseppe/local/bin/vcftools";

for FILE in ${PDATA}/*.vcf;
        do ID=`basename ${FILE} ".vcf"`;

        SCRIPT=vcftools_filter_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate' >>${PDATA}/${SCRIPT};
        echo "${PCODE} --vcf ${FILE} --remove-indels --maf ${MAF} --max-maf ${MAXMAF} --hwe ${HWE} --recode --recode-INFO-all --out ${PDATA}/${ID}_noindels_maf${MAF}_hwe${HWE}" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/vcftools_filter_${ID}.err -o ${PDATA}/vcftools_filter_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done

