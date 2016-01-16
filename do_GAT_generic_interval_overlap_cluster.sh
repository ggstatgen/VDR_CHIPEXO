#!/bin/bash
#script to run the GAT on one peak file at a time
#you can feed GAT with a segment bed where you cat'ed all the peakfiles, however using multiple processor is presumably faster than using multiple cores

#this should be used for a generic intersection of two intervals (eg peaks and DNASE)

#play around with the annotations and isochores etc

#also choose the KIND OF OVERLAP:
#Counters describe the measure of association that is tested. Counters are selected with the command line option --counter. Available counters are:
#
#    nucleotide-overlap: number of bases overlapping [default]
#    segment-overlap: number of intervals intervals in the segments of interest overlapping annotations. A single base-pair overlap is sufficient.
#    segment-mid-overlap: number of intervals in the segments of interest overlapping at their midpoint annotations.
#    annotations-overlap: number of intervals in the annotations overlapping segments of interest. A single base-pair overlap is sufficient.
#    segment-mid-overlap: number of intervals in the annotations overlapping at their midpoint segments of interest


if [ ! $# == 6 ]; then
        echo "Usage: `basename $0` <SEGMENT_PATH> <ANNOTATION> <MAPPABILITY> <NUM_SAMPLES> <OVERLAP> <DESC>"
        echo "<SEGMENT_PATH> path to directory containing the bed file(s)"
	echo "<ANNOTATION> full path to annotation bed file"
	echo "<MAPPABILITY> one of [36bp|40bp]"
	echo "<NUM_SAMPLES> number of samples to use for the randomisation (eg 50000)"
        echo "<OVERLAP> measure of association [nucleotide-overlap|segment-overlap]"
	echo "<DESC> short 1word description of the annotation used (eg wgEncodeRegDnaseClustered) - used for naming the out dir"
	echo "NOTE: workspace will be contigs_ungapped.bed.gz. Change if needed."
        exit
fi

PDATA=$1;
PANNOTATION=$2;
MAPPABILITY=$3;
PRAND=$4;
OVERLAP=$5;
DESC=$6;

#data
PISOCHORE="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/d_GAT/hg19.fasta.isochore.bed.gz";
#PISOCHORE="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/d_GAT/b37.fasta.isochore.bed.gz";

#generic workspace (B1 paper)
PWSP="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/d_GAT/contigs_ungapped.bed.gz";
#peak overlap workspace 2 (B2 paper)
#PWSP="/net/isi-scratch/giuseppe/CHIPSEQ_VDR/PAPER_PEAKS/vdr_chipseq_peaks_hg19_GAT_slop5kb.bed.gz";
#peak overlap workspace 3 (B3 paper)
#PWSP="/net/isi-scratch/giuseppe/indexes/Hsap/Ensembl73/ensembl73_known_protein_coding_genes_UCSC_sorted_plus5kbTSS_minus1kbTES_GATws.bed.gz";


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

#PMAPPABILITY="/net/isi-mirror/ucsc/hg19/mappability/40bp/hg19_40mer.bed.gz";
#PMAPPABILITY="/net/isi-mirror/ucsc/hg19/mappability/36bp/hg19_36mer.bed.gz";

PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";


POUT=${PDATA}/GAT_${DESC}_r${PRAND};
mkdir ${POUT};

for FILE in ${PDATA}/*.bed;
	do ID=`basename ${FILE} ".bed"`;

        SCRIPT=GAT_${DESC}_vs_${ID}_${PRAND}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
        echo 'source activate' >> ${PDATA}/${SCRIPT};

	#change required options here:
	#REMOVED               --num-threads=${THREADS}       \\
	echo "${PCODE}/gat-run.py     \\
              --ignore-segment-tracks        \\
              --segments=${FILE}             \\
              --counter=${OVERLAP}           \\
              --annotations=${PANNOTATION}   \\
              --workspace=${PWSP}            \\
              --workspace=${PMAPPABILITY}    \\
              --isochore-file=${PISOCHORE}   \\
              --output-counts-pattern=${POUT}/${ID}.%s.oc.tsv.gz \\
              --num-samples=${PRAND}         \\
              --log=${POUT}/GAT_${DESC}_vs_${ID}.log" >> ${PDATA}/${SCRIPT};
        nice -5 qsub  -v "BASH_ENV=~/.bashrc" -e ${POUT}/GAT_${DESC}_vs_${ID}.err -o ${POUT}/GAT_${DESC}_vs_${ID}.tsv -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
