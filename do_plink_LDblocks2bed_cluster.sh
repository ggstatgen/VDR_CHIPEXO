#!/bin/bash

#converts the output of plink --blocks to a bed file

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "<PATH> data path for the plink blocks.det files"
        exit
fi

PDATA=$1;
PCODE="/net/isi-backup/giuseppe/scripts/P_Various/do_LD_PLINK_blockfile2bed.pl";

for FILE in ${PDATA}/*.blocks.det;
        do ID=`basename ${FILE} ".blocks.det"`;

        SCRIPT=plink_LDblock2bed_${ID}.sh;
        
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate'  >>${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};
	echo "perl ${PCODE} -i=${FILE}" >>${PDATA}/${SCRIPT};

	nice -5 qsub -cwd -e plink_LDblock2bed_${ID}.err -o ${ID}_intervals.bed -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
