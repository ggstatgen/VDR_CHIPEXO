#!/bin/bash

#uses the c code kindly given to me by Daniel Taliun (see mail 2/12/2014)
#to try and derive genome wide LD blocks using an r^2 definition.
#it is based on the program LD explorer found here
#http://www.eurac.edu/en/research/health/biomed/services/Pages/LDExplorer.aspx

#Usage example:
#mig_rsq --vcf <vcf file with single chromosome> [--region <start bp> <end_bp>] --maf <maf threshold> --rsq <low LD> <high LD> --frac <fraction of low LD and high LD> --window 25 --blocks_out <output_file_1> <output_file_2>

#After --vcf option specify input VCF file.

#If you are interested in sub-region, then use optional --region command to specify start and end positions in base pairs. After --maf option specify minimal minor allele frequency. By default is 0. After --rsq specify two r^2 thresholds for lowLD and highLD. If r^2 < threshold_1 then LD between SNPs is weak, if r^2 >= threshold_2 then LD between SNPs is strong. All other situations are considered non-informative.

#After --frac specify a threshold for the fraction srtongLD/(strongLD + weakLD). I.e. the region is a block if at least frac of SNP pairs among all informative have r^2 >= threshold_2. E.g. if you set --rsq 0.5 0.5 --frac 0.95, then a region is a block if >0.95 of SNP pairs have r^2 >= 0.5.

#After --blocks_out specify two output files. Each output file contains a list of blocks. output_file_2 contains a list of all blocks (can be very long). output_file_1 contains a list of non-overlapping blocks selected from output_file_2 using greedy approach to optimize chromosome coverage and maximize blocks length.


#but its experimental and might not work well. Still testing

#trying the following test run
#/net/isi-scratch/giuseppe/tools/MIG_Rsq/src/mig_rsq --vcf /net/isi-mirror/1000_Genomes/POPULATION_PANELS/1kg_VCF/CEU/chr19.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf --maf 0.002 --rsq 0.5 0.5 --frac 0.95 --window 25 --blocks_out test1 test2 &> rsquare_test.out &


#MAF=0.002;
#FRAC=0.95;
WINDOW=25;

if [ ! $# == 4 ]; then
        echo "Usage: `basename $0` <PATH> <LD_THRS> <FRAC> <MAF>"
        echo "<PATH> data path for the vcf (NOT .gz) files to compute LD blocks from"
	echo "<LD_THRS> r^2 LD threshold"
	echo "<FRAC> the region is a block if at least FRAC of SNP pairs have r^2 LD_THRES"
	echo "<MAF> minimum MAF  to consider (eg 0.01)"
        exit
fi

PDATA=$1;
LD=$2;
FRAC=$3;
MAF=$4;
PCODE="/net/isi-scratch/giuseppe/tools/MIG_Rsq/src";

for FILE in ${PDATA}/*.vcf;
        do ID=`basename ${FILE} ".vcf"`;

        SCRIPT=CEU_MIGLD_r${LD}_MAF${MAF}_FRAC${FRAC}_${ID}.sh;
        
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate'  >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "${PCODE}/mig_rsq --vcf ${FILE} --maf ${MAF} --rsq ${LD} ${LD} --frac ${FRAC} --window ${WINDOW} --blocks_out ${PDATA}/MIG_${ID}_${LD}.txt ${PDATA}/MIG_${ID}_${LD}.txt.gz" >>${PDATA}/${SCRIPT};

	nice -5 qsub -cwd -e MIGLD_r${LD}_MAF${MAF}_FRAC${FRAC}_${ID}.err -o MIGLD_r${LD}_MAF${MAF}_FRAC${FRAC}_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT}; 
done
