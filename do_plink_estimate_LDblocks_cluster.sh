#!/bin/bash

#uses plink to estimate ld blocks given plink bed data (see do_vcf2plink...)

#commands taken from
#Getting LD blocks from 1000g using PLINK and a pipeline found here
#https://www.biostars.org/p/2909/#75824 ( answer 4)

#vcftools --gzvcf genotypes.subset.vcf.gz --plink-tped --out plinkformat
#plink --tfile plinkformat --make-bed --out plinkBEDformat --noweb --maf 0.002 --hwe 0.001
#mv plinkBEDformat.fam plinkBEDformat.fam.tmp
#cat plinkBEDformat.fam.tmp | sed 's/-9/1/g' > plinkBEDformat.fam

#time plink --bfile plinkBEDformat  --r2 inter-chr --ld-window-r2 $(RSQUARE) --out plinkoutput --noweb --threads 10 
#or maybe
#plink --bfile mydata --blocks 

#see here
#http://pngu.mgh.harvard.edu/~purcell/plink/ld.shtml

#and here for plink 1.9
#https://www.cog-genomics.org/plink2/ld

#documentation for blocks
#PLINK v1.90b2m 64-bit (15 Oct 2014)        https://www.cog-genomics.org/plink2
#(C) 2005-2014 Shaun Purcell, Christopher Chang   GNU General Public License v3
#Logging to plink.log.
#Error: No input dataset.
#For more information, try 'plink --help [flag name]' or 'plink --help | more'.
#(ve)[giuseppe@fgu217 d_VCF]$ plink --help blocks
#PLINK v1.90b2m 64-bit (15 Oct 2014)        https://www.cog-genomics.org/plink2
#(C) 2005-2014 Shaun Purcell, Christopher Chang   GNU General Public License v3
#
#--blocks <no-pheno-req> <no-small-max-span>
#  Estimate haplotype blocks, via Haploview's interpretation of the block
#  definition suggested by Gabriel S et al. (2002) The Structure of Haplotype
#  Blocks in the Human Genome.
#  * Normally, individuals with missing phenotypes are not considered by this
#    computation; the 'no-pheno-req' modifier lifts this restriction.
#  * Normally, size-2 blocks may not span more than 20kb, and size-3 blocks
#    are limited to 30kb.  The 'no-small-max-span' modifier removes these
#    limits.
#  The .blocks file is valid input for PLINK 1.07's --hap command.  However,
#  the --hap... family of flags has not been reimplemented in PLINK 1.9 due to
#  poor phasing accuracy relative to other software; for now, we recommend
#  using BEAGLE instead of PLINK for case/control haplotype association
#  analysis.  (You can use '--recode beagle' to export data to BEAGLE 3.3.)
#  We apologize for the inconvenience, and plan to develop variants of the
#  --hap... flags which handle pre-phased data effectively.
#
#--blocks-max-kb [kbs]      : Set --blocks maximum haploblock span (def. 200).
#--blocks-min-maf [cutoff]  : Adjust --blocks MAF minimum (default 0.05).
#--blocks-strong-lowci [x]  : Set --blocks 'strong LD' CI thresholds (defaults
#--blocks-strong-highci [x]   0.70 and 0.98).
#--blocks-recomb-highci [x] : Set 'recombination' CI threshold (default 0.90).
#--blocks-inform-frac [x]   : Force haploblock [strong LD pairs]:[total
#                             informative pairs] ratios to be larger than this
#                             value (default 0.95).


#more documentation: https://www.cog-genomics.org/plink2/ld
#based on what indicated here, I decided to use the option --blocks-strong-lowci 0.7005


if [ ! $# == 5 ]; then
        echo "Usage: `basename $0` <PATH> <LOWCI> <HICI> <RECOMBCI> <INFORMFRAC>"
        echo "<PATH> data path for the plink bed,bim,fam files to convert"
	echo "<LOWCI> --blocks strong LD CI low threshold (default 0.7005)"
	echo "<HICI> --blocks strong LD CI high threshold (default 0.98)"
	echo "<RECOMBCI> recombination CI threshold (default 0.90)"
	echo "<INFORMFRAC> (strong LD pairs)/(total informative pairs) > than this value (default 0.95)"
        exit
fi

PDATA=$1;
LOWCI=$2;
HICI=$3;
RECOMBCI=$4;
INFORMFRAC=$5;

PCODE="/home/giuseppe/local/bin";
BLOCKS_MAX_KB=500; #for coherency with SNAP (which uses a 500KB distance limit)

for FILE in ${PDATA}/*.bed;
        do ID=`basename ${FILE} ".bed"`;
        #ID=`echo ${FILE} | egrep -o "chr[123456789XYM]*"`;

        SCRIPT=plinkLD_${ID}_lowci_${LOWCI}_hici_${HICI}_recombci_${RECOMBCI}_informfrac_${INFORMFRAC}.sh;
        
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate'  >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	#echo "${PCODE}/plink --bfile ${PDATA}/${ID} --blocks  no-pheno-req --blocks-strong-lowci 0.7005 --blocks-max-kb ${BLOCKS_MAX_KB} --out ${PDATA}/plink_${ID} --threads 15" >>${PDATA}/${SCRIPT};
	#echo "${PCODE}/plink --bfile ${PDATA}/${ID} --blocks  no-pheno-req --blocks-strong-lowci 0.6005 --blocks-strong-highci 0.8305 --blocks-recomb-highci 0.7 --blocks-max-kb ${BLOCKS_MAX_KB} --out ${PDATA}/plink_${ID} --threads 15" >>${PDATA}/${SCRIPT};

	echo "${PCODE}/plink --bfile ${PDATA}/${ID} --blocks  no-pheno-req --blocks-strong-lowci ${LOWCI} --blocks-strong-highci ${HICI} --blocks-recomb-highci ${RECOMBCI} --blocks-inform-frac ${INFORMFRAC} --blocks-max-kb ${BLOCKS_MAX_KB} --out ${PDATA}/plink_${ID}_lowci_${LOWCI}_hici_${HICI}_recombci_${RECOMBCI}_informfrac_${INFORMFRAC} --threads 10" >>${PDATA}/${SCRIPT};

	nice -5 qsub -cwd -e plinkLD_${ID}_lowci_${LOWCI}_hici_${HICI}_recombci_${RECOMBCI}_informfrac_${INFORMFRAC}.err -o plinkLD_${ID}_lowci_${LOWCI}_hici_${HICI}_recombci_${RECOMBCI}_informfrac_${INFORMFRAC}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
