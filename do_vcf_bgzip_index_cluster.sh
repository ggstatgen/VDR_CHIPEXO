#!/bin/bash

#compress using bgzip
#index using tabix


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
        echo "You must specify the data path for the vcf files (e.g. /home/me/files)"
        exit
fi

PDATA=$1;
PCODE="/home/giuseppe/local/bin";

for FILE in ${PDATA}/*.vcf;
        do ID=`basename ${FILE} ".vcf"`;

        SCRIPT=vcf_bgzip_tabix_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo "${PCODE}/bgzip ${FILE}" >>${PDATA}/${SCRIPT};
	echo "${PCODE}/tabix -p vcf ${PDATA}/${ID}.vcf.gz" >>${PDATA}/${SCRIPT};

	nice -5 qsub -cwd -e vcf_bgzip_tabix_${ID}.err -o vcf_bgzip_tabix_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
