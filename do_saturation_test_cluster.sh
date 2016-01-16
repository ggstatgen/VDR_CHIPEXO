#!/bin/bash
#this calls the perl script saturation_Test_mt.pl
#and shares the workload on the cluster

#USAGE: do_saturation_test.pl -input=<INFILE> -bin=<BIN>
#<INFILE> input fastq file
#<BIN>: number of mapped reads in bin (e.g. 100000)



if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH_DATA> <BIN>"
        echo "PATH_DATA - directory containing the input .fastq files"
	echo "BIN - size of the bin in number of mapped reads (e.g. 1000000, sample every 1million reads)"
        exit
fi

PDATA=$1;
PCODE="/net/isi-backup/giuseppe/scripts";

for FILE in ${PDATA}/*.fastq;
	do ID=`basename ${FILE} ".fastq"`;
        #do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

	SCRIPT=saturation_${ID}.sh;

	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate' >>${PDATA}/${SCRIPT};
	echo "perl ${PCODE}/encode_analyses.pl -i=${FILE} -bin=$2" >>${PDATA}/${SCRIPT};
        nice -5 qsub -v "BASH_ENV=~/.bashrc" -e ${PDATA}/saturation_${ID}.err -o ${PDATA}/saturation_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
