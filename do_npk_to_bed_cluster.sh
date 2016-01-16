#!/bin/bash
#31/10 MACS2 does not output beds anymore.
#This converts its narrowPeak files to bed

#USAGE: perl do_npk_to_bed.pl -i=<FILE> -b=<N_BEST>
#<FILE> input narrowPeak file
#<N_BEST> number of best peaks (by p-value) to return (optional, default = all)

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <N_BEST>"
        echo "<PATH> path for the .narrowPeak files"
	echo "<N_BEST> number of best peaks (ordered by p-val) to retrieve. See to 0 to retrieve all peaks"
        exit
fi

PDATA=$1;
PBEST=$2;
PCODE="/net/isi-backup/giuseppe/scripts";

for FILE in ${PDATA}/*.narrowPeak;
	do ID=`basename ${FILE} ".narrowPeak"`;

        SCRIPT=npk2bed_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	if [ ${PBEST} -eq 0 ]
	then
		echo "perl $PCODE/do_npk_to_bed.pl -i=${FILE}" >> ${PDATA}/${SCRIPT};
	else
		echo "perl $PCODE/do_npk_to_bed.pl -i=${FILE} -b=${PBEST}" >> ${PDATA}/${SCRIPT};
	fi

        nice -5 qsub -e ${PDATA}/npk2bed_${ID}.err -o ${PDATA}/${ID}.bed -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
