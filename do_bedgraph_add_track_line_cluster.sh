#!/bin/bash

#30/5/2014
#add a track line to the top of each bedgraph

if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <EXT>"
        echo "PATH - path for the begraph files"
	echo "EXT - file extension for the files"
        exit
fi

PDATA=$1;
EXT=$2;

for FILE in ${PDATA}/*.${EXT};
        do 
	BASENAME=`basename ${FILE} ".${EXT}"`;
        ID=`echo ${FILE} | egrep -o "NA[0-9]*"`;
	SCRIPT=bg_addtrack_${ID}.sh;

	#LINE="track type=bedGraph name=\\\"${ID}\\\" description=\\\"WIGGLER-normalised track for sample ${ID}\\\""
	LINE="track type=bedGraph name=\\\"${ID}\\\" description=\\\"WIGGLER-normalised track for sample ${ID}\\\" autoScale=off windowingFunction=mean+whiskers  maxHeightPixels=100:100:5 smoothingWindow=5"

	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
	#echo "awk '{print \$1, \"\t\", \$2, \"\t\", \$3}' ${FILE}" >>${PDATA}/${SCRIPT};
	echo "awk 'BEGIN{print \"${LINE}\"}{print}' ${FILE} > ${PDATA}/${BASENAME}_h.${EXT}" >>${PDATA}/${SCRIPT};

	#add a bedgraph header
	#awk 'BEGIN{print "new line"}{print}' infile > outfile
	nice -5 qsub -e ${PDATA}/bg_addtrack_${ID}.err -o ${PDATA}/bg_addtrack_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT};
done

#add this stuff to the WIGGLER script so it does it automatically
