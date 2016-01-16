#!/bin/bash
#31/10 MACS2 does not output beds anymore.
#This converts its narrowPeak files to bed

#USAGE: perl do_npk_to_bed.pl -i=<FILE> -b=<N_BEST>
#<FILE> input narrowPeak file
#<N_BEST> number of best peaks (by p-value) to return (optional, default = all)

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <EXT>"
        echo "<PATH> path for the .narrowPeak files"
	echo "<EXT> file extension"
	echo "collects chr, start, stop, name, pval"
        exit
fi

PDATA=$1;
EXT=$2;

for FILE in ${PDATA}/*.${EXT};
	do ID=`basename ${FILE} ".${EXT}"`;

        SCRIPT=npk2bed_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	
	echo "cat ${FILE} | cut -f 1,2,3,4,8" >> ${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/npk2bed_${ID}.err -o ${PDATA}/${ID}.bed -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
