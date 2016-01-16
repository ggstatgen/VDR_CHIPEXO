#!/bin/bash

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <EXTENSION>"
        echo "<PATH> Data path for the interval files (e.g. /home/me/files)"
	echo "<EXTENSION> file extension for the interval files (bed, npk, etc)"
        exit
fi

PDATA=$1;
EXT=$2;
#location of of the main script
PCODE="/net/isi-backup/giuseppe/scripts";
PERL="/net/isi-cgat/ifs/apps/apps/perl-5.16.1/bin/perl";

for FILE in ${PDATA}/*.${EXT};
	do ID=`basename ${FILE} ".${EXT}"`;

        SCRIPT=pklength_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
#	echo "source activate" >> ${PDATA}/${SCRIPT};
	echo "${PERL} ${PCODE}/do_bed_get_peak_sizes.pl -i=${FILE} > ${PDATA}/pklenght_${ID}.hist" >>${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/pklength_${ID}.err -o ${PDATA}/pklength_${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
