#!/bin/bash
#script to perform library complexity estimation using Picard on the FGU clusters
#Giuseppe Gallone 14/2/2013

#ONLY WORKS WITH PAIRED END SO USELESS with VDR

#picard docs
# Attempts to estimate library complexity from sequence of read pairs alone.

#MIN_IDENTICAL_BASES=Integer	
#The minimum number of bases at the starts of reads that must be identical for reads to be grouped together for duplicate detection. 
#In effect total_reads / 4^max_id_bases reads will be compared at a time, so lower numbers will produce more accurate results but consume exponentially more memory and CPU. 
#Default value: 5. This option can be set to 'null' to clear the default value.

#MAX_DIFF_RATE=Double	
#The maximum rate of differences between two reads to call them identical. Default value: 0.03. This option can be set to 'null' to clear the default value.

#MIN_MEAN_QUALITY=Integer	
#The minimum mean quality of the bases in a read pair for the read to be analyzed. Reads with lower average quality are filtered out and not considered in any calculations. 
#Default value: 20. This option can be set to 'null' to clear the default value.

#MAX_GROUP_RATIO=Integer
#Do not process self-similar groups that are this many times over the mean expected group size. 
#I.e. if the input contains 10m read pairs and MIN_IDENTICAL_BASES is set to 5, then the mean expected group size would be approximately 10 reads. 
#Default value: 500. This option can be set to 'null' to clear the default value.


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input .bam files"
	exit
fi

PDATA=$1;

#------------------------
#hard coded paths
#------------------------
PCODE_PICARD="/net/isi-scratch/giuseppe/tools/picard-tools-1.117";
TMP="/tmp";


for FILE in ${PDATA}/*.bam;
        do
        FILE_ID=`basename ${FILE} ".bam"`;

        SCRIPT=picard_elc_${FILE_ID}.sh;
        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        echo "java -Xmx4g -Djava.io.tmpdir=${TMP} \\
              -jar ${PCODE_PICARD}/EstimateLibraryComplexity.jar \\
              INPUT=${FILE} \\
              OUTPUT=${PDATA}/${FILE_ID}.metrics" >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	nice -5 qsub -e ${PDATA}/picard_elc_${FILE_ID}.stderr -o ${PDATA}/picard_elc_${FILE_ID}.stdout -v "BASH_ENV=~/.bashrc" -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done
