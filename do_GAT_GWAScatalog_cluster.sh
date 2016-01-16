#!/bin/bash

#script to run the GAT on one peak file at a time
#you can feed GAT with a segment bed where you cat'ed all the peakfiles, however using multiple processor is presumably faster than using multiple cores

#play around with the annotations and isochores etc

#also choose the KIND OF OVERLAP:
#Counters describe the measure of association that is tested. Counters are selected with the command line option --counter. Available counters are:
#
#    nucleotide-overlap: number of bases overlapping [default]
#    segment-overlap: number of intervals intervals in the segments of interest overlapping annotations. A single base-pair overlap is sufficient.
#    segment-mid-overlap: number of intervals in the segments of interest overlapping at their midpoint annotations.
#    annotations-overlap: number of intervals in the annotations overlapping segments of interest. A single base-pair overlap is sufficient.
#    segment-mid-overlap: number of intervals in the annotations overlapping at their midpoint segments of interest

#CHAT WITH ANDREAS
#IMPORTANT OPTION
#  --conditional=CONDITIONAL
#                        conditional workspace creation [default=unconditional]
#                        cooccurance - compute enrichment only within workspace
#                        segments that contain both segments  and annotations
#                        annotation-centered - workspace centered around
#                        annotations. See --conditional-extension segment-
#                        centered - workspace centered around segments. See
#                        --conditional-extension
#  --conditional-extension=CONDITIONAL_EXTENSION
#                        if workspace is created conditional, extend by this
#                        amount (in bp) [default=none].

#you need to run with --conditional=annotation-centered
#Imagine some annotation which cover a significant part of the genome and your peaks all fall in it. Then of course you will see some enrichment.
#you should restrict your workspace to only the annotation region

#OUTPUT COUNTS
#for example useful when you want to compare the enrichments in two GAT runs
#--output-counts-pattern=CD4.%s.overlap.counts.tsv.gz


if [ ! $# == 6 ]; then
        echo "Usage: `basename $0` <SEGMENT_PATH> <ANNOTATION> <MAPPABILITY> <NUM_SAMPLES> <OVERLAP> <DESC>"
        echo "<SEGMENT_PATH> path to directory containing the bed file(s)"
        echo "<ANNOTATION> full path to annotation bed file"
	echo "<MAPPABILITY> one of [36bp|40bp]"
        echo "<NUM_SAMPLES> number of samples to use for the randomisation (eg 10000)"
        echo "<OVERLAP> measure of association [nucleotide-overlap|segment-overlap]"
        echo "<DESC> short 1word description of the annotation used (eg 8snp_150kb) - used for naming the out dir"
        exit
fi

PDATA=$1;
PANNOTATION=$2;
MAPPABILITY=$3;
PRAND=$4;
OVERLAP=$5;
DESC=$6;

#data
#PANNOTATION="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/d_GAT/annotations_geneset.bed.gz";
#PANNOTATION="/net/isi-scratch/giuseppe/indexes/GWAS/DistiLDsnps_segments.bed";
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

POUT=${PDATA}/GAT_GWAScatalog_r${PRAND}_a${DESC};
mkdir ${POUT};

for FILE in ${PDATA}/*.bed;
	do ID=`basename ${FILE} ".bed"`;

        SCRIPT=GAT_GWAScat_${ID}_${PRAND}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'source activate' >> ${PDATA}/${SCRIPT};
	echo '' >>${PDATA}/${SCRIPT};


	#change required options here:
	#REMOVED               --num-threads=${THREADS}       \\
	#--conditional=annotation-centered \\
        #--conditional-extension=100
	echo "${PCODE}/gat-run.py     \\
              --ignore-segment-tracks        \\
              --segments=${FILE}             \\
              --counter=${OVERLAP}           \\
              --annotations=${PANNOTATION}   \\
              --isochore-file=${PISOCHORE}   \\
              --workspace=${PWSP}            \\
              --workspace=${PMAPPABILITY}   \\
              --output-counts-pattern=${POUT}/${ID}.%s.oc.tsv.gz \\
              --num-samples=${PRAND}         \\
              --log=${POUT}/GAT_GWAScat_vs_${ID}.log" >> ${PDATA}/${SCRIPT};
        nice -5 qsub  -v "BASH_ENV=~/.bashrc" -e ${POUT}/GAT_GWAScat_vs_${ID}.err -o ${POUT}/GAT_GWAScat_vs_${ID}.tsv -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
