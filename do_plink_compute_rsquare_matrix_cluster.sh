#!/bin/bash

#uses plink to produce r^2 matrices 
#from these you should be able to obtain  r^2 intervals as you please
#this follows a mail from the author, who says I cannot use blocks to estimate this
#No, --blocks is always based on D' confidence intervals.  However, you can write a custom script that postprocesses the output of --r2.


#more documentation: https://www.cog-genomics.org/plink2/ld
#based on what indicated here, I decided to use the option --blocks-strong-lowci 0.7005

#--r2 <square | square0 | triangle | inter-chr> <gz | bin> <single-prec>
#     <spaces> <in-phase> <dprime> <with-freqs> <yes-really>
#  LD statistic reports.  --r yields raw inter-variant correlations, while
#  --r2 reports their squares.  You can request results for all pairs in
#  matrix format (if you specify 'bin' or one of the shape modifiers), all
#  pairs in table format ('inter-chr'), or a limited window in table format
#  (default).
#  * The 'gz' modifier causes the output text file to be gzipped.
#  * 'bin' causes the output matrix to be written in binary format.  The
#    matrix is square if no shape is explicitly specified.
#  * In combination with 'bin', 'single-prec' causes single-precision instead
#    of double-precision numbers to be written.
#  * By default, text matrices are tab-delimited; 'spaces' switches this.
#  * 'in-phase' adds a column with in-phase allele pairs to table-formatted
#    reports.  (This cannot be used with very long allele codes.)
#  * 'dprime' adds Lewontin's D-prime statistic to table-formatted reports,
#    and forces both r/r^2 and D-prime to be based on the maximum likelihood
#    solution to the cubic equation discussed in Gaunt T, Rodriguez S, Day I
#    (2007) Cubic exact solutions for the estimation of pairwise haplotype
#    frequencies.
#  * 'with-freqs' adds MAF columns to table-formatted reports.
#  * Since the resulting file can easily be huge, you're required to add the
#    'yes-really' modifier when requesting an unfiltered, non-distributed all
#    pairs computation on more than 400k variants.
#  * These computations can be subdivided with --parallel (even when the
#    'square' modifier is active).

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "<PATH> data path for the plink bed,bim,fam files to convert"
        exit
fi

PDATA=$1;
PCODE="/home/giuseppe/local/bin";
BLOCKS_MAX_KB=500; #for coherency with SNAP (which uses a 500KB distance limit)

for FILE in ${PDATA}/*.bed;
        do ID=`basename ${FILE} ".bed"`;
        #ID=`echo ${FILE} | egrep -o "chr[123456789XYM]*"`;

        SCRIPT=plink_compute_rsquares_${ID}.sh;
        
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate'  >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "${PCODE}/plink --bfile ${PDATA}/${ID} --r2 gz dprime yes-really --out ${PDATA}/plink_${ID} --threads 15"  >>${PDATA}/${SCRIPT};
	#echo "${PCODE}/plink --bfile ${PDATA}/${ID} --r2 inter-chr gz dprime yes-really --out ${PDATA}/plink_${ID} --threads 15"  >>${PDATA}/${SCRIPT};
	#echo "${PCODE}/plink --bfile ${PDATA}/${ID} --blocks  no-pheno-req --blocks-strong-lowci 0.7005 --blocks-max-kb ${BLOCKS_MAX_KB} --out ${PDATA}/plink_${ID} --threads 15" >>${PDATA}/${SCRIPT};
	nice -5 qsub -cwd -e plink_compute_rsquares_${ID}.err -o plink_compute_rsquares_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
