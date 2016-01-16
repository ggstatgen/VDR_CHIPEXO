#!/bin/bash

#24/04
#Script to run MACS2+IDR pseudoreplicate peak calling on the cluster
# The script returns a list of MACS2-called peaks thresholded via IDR
# IDR is performed by splitting the file in two random halfs (pseudo replicates) as suggested by Anshul Kundake 
# this is due to the absence of controls or real pseudoreps
#The output bedfile should be a set of strong MACS2 peaks which passed IDR control during pseudorep comparison

#IDR threshold and macs 2 parameters are hardcoded in the perl script. Change that

if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "You must specify the data path for the bam files (e.g. /home/me/files)"
        exit
fi

PDATA=$1;
PCODE="/net/isi-backup/giuseppe/scripts";

for FILE in ${PDATA}/*.bam;
	do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

        SCRIPT=${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	
	echo "perl ${PCODE}/do_macs2_idr_pseudoreps.pl -i=${FILE}" >>${PDATA}/${SCRIPT};
        nice -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
