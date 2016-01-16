#!/bin/bash

#script to run the GAT against GENOMIC REGIONS (introns,exons, etc) on one peak file at a time
#you can feed GAT with a segment bed where you cat'ed all the peakfiles, however using multiple processor is presumably faster than using multiple cores

#when the segments are too large, GAT complains
#play with the following parameters

#  --bucket-size=BUCKET_SIZE
#                        size of a bin for histogram of segment lengths
#                        [default=1]
#  --nbuckets=NBUCKETS   number of bins for histogram of segment lengths
#                        [default=100000]



#play around with the annotations and isochores etc


if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <NUM_SAMPLES> <MAPPABILITY>"
        echo "<PATH> data path for the bed segment file"
	echo "<NUM_SAMPLES> number of samples to use for the randomisation (10000)"
        echo "<MAPPABILITY> one of [36bp|40bp]"
        exit
fi

PDATA=$1;
PRAND=$2;
MAPPABILITY=$3;

#data
#PANNOTATION="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/d_GAT/annotations_geneset.bed.gz";
PANNOTATION="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/d_GAT/annotation_hg19_ens72.bed.gz";
PISOCHORE="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/d_GAT/hg19.fasta.isochore.bed.gz";
PWSP="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/d_GAT/contigs_ungapped.bed.gz"

if [ "$MAPPABILITY" = "36bp" ]
then
        PMAPPABILITY="/net/isi-mirror/ucsc/hg19/mappability/36bp/hg19_36mer.bed.gz";
elif [ "$MAPPABILITY" = "40bp" ]
then
        PMAPPABILITY="/net/isi-mirror/ucsc/hg19/mappability/40bp/hg19_40mer.bed.gz";
else
        echo "ERROR: mappability field not recognised. Exiting.";
        exit 1;
fi

#PMAPPABILITY="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/d_GAT/mapability_36.filtered.bed.gz";
#PMAPPABILITY="/net/isi-mirror/ucsc/hg19/mappability/40bp/hg19_40mer.bed.gz";


PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";
POUT=${PDATA}/GAT_GRA_r${PRAND};
mkdir ${POUT};

for FILE in ${PDATA}/*.bed;
	do ID=`basename ${FILE} ".bed"`;

        SCRIPT=GAT_gra_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate' >>${PDATA}/${SCRIPT};
	#change required options here:
	#removed -t ${PTHREADS}               \\
	echo "${PCODE}/gat-run.py          \\
              --ignore-segment-tracks      \\
              --segments=${FILE}           \\
              --isochore-file=${PISOCHORE} \\
              --workspace=${PMAPPABILITY}  \\
              --annotations=${PANNOTATION} \\
              --workspace=${PWSP}          \\
              --num-samples=${PRAND}       \\
              --output-counts-pattern=${POUT}/${ID}.%s.oc.tsv.gz \\
              --bucket-size=5              \\
              --nbuckets=2000000            \\
              --log=${POUT}/GAT_${ID}.log" >> ${PDATA}/${SCRIPT};
        
        nice -5 qsub -e ${POUT}/GAT_gra_${ID}.err  -o ${POUT}/GAT_gra_${ID}.tsv -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
