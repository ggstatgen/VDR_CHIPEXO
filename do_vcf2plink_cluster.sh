#!/bin/bash

#uses vcf-tools to convert from vcf.gz to plink tfamm then to plink bed

#commands taken from
#Getting LD blocks from 1000g using PLINK and a pipeline found here
#https://www.biostars.org/p/2909/#75824 ( answer 4)

#vcftools --gzvcf genotypes.subset.vcf.gz --plink-tped --out plinkformat
#plink --tfile plinkformat --make-bed --out plinkBEDformat --noweb --maf 0.002 --hwe 0.001
#mv plinkBEDformat.fam plinkBEDformat.fam.tmp
#cat plinkBEDformat.fam.tmp | sed 's/-9/1/g' > plinkBEDformat.fam


if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <MAF>"
        echo "<PATH> data path for the vcf.gz files to convert"
	echo "<MAF> minimum maf to retain (eg 0.002)"
        exit
fi

PDATA=$1;
MAF=$2;
PCODE="/home/giuseppe/local/bin";

for FILE in ${PDATA}/*.vcf.gz;
        do ID=`basename ${FILE} ".vcf.gz"`;
        #do ID=`echo ${FILE} | egrep -o "chr[123456789XYM]*"`;

        SCRIPT=vcf2plink_${ID}_MAF${MAF}.sh;
        
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate'  >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "${PCODE}/vcftools --gzvcf ${FILE} --plink-tped --out ${PDATA}/${ID}" >>${PDATA}/${SCRIPT};
	echo "${PCODE}/plink --tfile ${PDATA}/${ID} --make-bed --out ${PDATA}/${ID} --noweb --maf ${MAF} --hwe 0.001" >>${PDATA}/${SCRIPT};
	echo "mv ${PDATA}/${ID}.fam ${PDATA}/${ID}.fam.tmp" >>${PDATA}/${SCRIPT};
	echo "cat ${PDATA}/${ID}.fam.tmp | sed 's/-9/1/g' > ${PDATA}/${ID}.fam" >>${PDATA}/${SCRIPT};

	nice -5 qsub -cwd -e vcf2plink_${ID}_MAF${MAF}.err -o vcf2plink_${ID}_MAF${MAF}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT}; 
done
