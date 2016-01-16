#!/bin/bash


#script to post-process the output of the MIGLD r^2 algorithms for finding r^2 based LD blocks
#this will get beds from the output with
#bed name: LD BLOCK ID
#bed score: number of markers in the BLOCK

#cat MIG_chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.AFR_0.9.txt | awk '{OFS="\t"; if (!/^#/){print "chr21",$6,$7,$1,$8}}'  | sort -k1,1V -k2,2g | uniq > MIG_chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.AFR_0.9.bed

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "<PATH> data path for the MIG .txt files"
        exit
fi

PDATA=$1;

#MIG_chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.AFR_0.9.txt

for FILE in ${PDATA}/MIG*.txt;
	do 
	BASENAME=`basename ${FILE} ".txt"`;
	ID=`expr "$FILE" : '.*_\(chr[X0-9]*\).'`;

        SCRIPT=MIGLD_to_bed_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "cat ${FILE} | awk '{OFS=\"\t\"; if (!/^#/ && !/^BLOCK_NAME/){print \"$ID\",\$6,\$7,\$1,\$8}}'  | sort -k1,1V -k2,2g | uniq > ${PDATA}/${BASENAME}.bed" >>${PDATA}/${SCRIPT};
	
        nice -5 qsub -e ${PDATA}/MIGLD_2_bed_${ID}.err -o ${PDATA}/MIGLD_2_bed_${ID}.out  -q short_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
