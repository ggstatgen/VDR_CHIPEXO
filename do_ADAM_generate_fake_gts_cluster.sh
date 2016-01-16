#!/bin/bash

#1/10/2015
#This processes some input vcf files from 1kg to remove the genotype info and replace it with "fake" genotype info to create a FULL ALTERNATE reference genome (using Alleleseq's vcf2diploid) and a FULL REFERENCE reference genome.
#Adam want the full alternate to run FIMO and test the presence of footprints/motifs


if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <INDELS>"
	echo "<INDELS> one of [YES|NO]; if no, only SNPs will be kept";
        exit
fi

PDATA=$1;
INDELS=$2;
PERL="/net/isi-cgat/ifs/apps/apps/perl-5.16.1/bin/perl";

PCODE="/net/isi-backup/giuseppe/scripts/P_Various/do_ADAM_generate_fake_gts.pl";

for FILE in ${PDATA}/chr*.vcf.gz;
        do
	ID=`basename ${FILE} ".vcf.gz"`;

	if [ "$INDELS" = "NO" ]
	then
	        SCRIPT=vcfscript_noindels_${ID}.sh;
	elif [ "$INDELS" = "YES" ]
	then
		SCRIPT=vcfscript_withindels_${ID}.sh
	else
	        echo "ERROR: INDELS field not recognised. Aborting.";
	        exit 1;
	fi

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate' >>${PDATA}/${SCRIPT};

	if [ "$INDELS" = "NO" ]
        then
		echo "$PERL ${PCODE} -i=${FILE} -snpsonly" >>${PDATA}/${SCRIPT};
        	nice -5 qsub -v "BASH_ENV=~/.bashrc" -e ${PDATA}/vcfscript_noindels_${ID}.log -o ${PDATA}/vcfscript_${ID}_noindels.vcf -q newnodes.q ${PDATA}/${SCRIPT};
	elif [ "$INDELS" = "YES" ]
        then
                echo "$PERL ${PCODE} -i=${FILE}" >>${PDATA}/${SCRIPT};
                nice -5 qsub  -v "BASH_ENV=~/.bashrc" -e ${PDATA}/vcfscript_withindels_${ID}.log -o ${PDATA}/vcfscript_${ID}_withindels.vcf -q newnodes.q ${PDATA}/${SCRIPT};                
        else
                echo "ERROR: INDELS field not recognised. Aborting.";
                exit 1;
        fi

        rm ${PDATA}/${SCRIPT};  
done
