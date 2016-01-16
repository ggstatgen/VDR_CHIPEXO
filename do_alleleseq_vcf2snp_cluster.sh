#!/bin/bash


#script to convert vcf files to SNP files required by the alleleseq pipeline
#http://info.gersteinlab.org/AlleleSeq

#Usage:
#     vcf2snp <vcf> > out.txt
#
#     -s default 1; only take in SNPs, remove all indels - set to 0 if indels need to be included
#     -p default 1; only take in SNPs with FILTER field PASS - if the field is non-pass for ALL snps, turn this off.
#     -c child ; sample ID in VCF file (mandatory)
#     -m mom ; if sample ID not given, assumed to be a non-trio conversion (no parents)
#     -d dad
#     -r default 0; to remove homozygous ref and alt SNPs, set to 1
#     -h help
#
#      Convert vcf format to file needed by alleleseq pipeline
#
#      output is:
#      chr     pos     ref_allele     Mat     Pat     Child   Phase
#
#      NOTE:
#      For the NONtrio option,
#      whenever the genotype is not 0|0 or 0/0 but like this 1:2.000:-5.00,-0.00, the script treat it as missing and skips it
#
#      Example:
#         zcat snp.vcf.gz | vcf2snp -p 0 -c NA12878 -m NA12892 -d NA12891 - > snp.call

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <SNPS_ONLY>"
        echo "<PATH> data path for the vcf.gz files"
	echo "<SNPS_ONLY> [1|0]. 1 = take only snps from vcf. 0 = take everything."
	echo "NOTE: this works with vcf containing ONLY the child, not trios"
        exit
fi

PDATA=$1;
SNPS_ONLY=$2;
#location of macs2 executable, change if needed
PCODE="/net/isi-scratch/giuseppe/tools/AlleleSeq_pipeline_v1.2";

for FILE in ${PDATA}/*.vcf.gz;
	do 
	#ID=`basename ${FILE} ".vcf.gz"`;
	ID=`echo ${FILE} | egrep -o "NA[0-9]*"`;

        SCRIPT=vcf2snp_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "source activate" >> ${PDATA}/${SCRIPT};

	echo "zcat ${FILE} | ${PCODE}/vcf2snp -s ${SNPS_ONLY} -p 0 -c ${ID} -"  >> ${PDATA}/${SCRIPT};
       
        nice -5 qsub -e ${PDATA}/vcf2nsp_${ID}.err -o ${PDATA}/vcf2snp_${ID}.snv -q newnodes.q ${PDATA}/${SCRIPT};
        #rm ${PDATA}/${SCRIPT};  
done
