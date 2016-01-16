#!/bin/bash

#This is to be used when the segments AND annotations are SNPs.
#Ultimately it is like doing a series of hypergeometric tests. Though GAT allows for background correction. 
#So if your segment SNPs are ASB snps, for example, they should be corrected because they will tend to be under peaks and therefore will be biased for GC-rich rregions?

#overlap should be segment-overlap. It should not make a difference though

#these are the possible kinds of overlap
#Counters describe the measure of association that is tested. Counters are selected with the command line option --counter. Available counters are:
#
#    nucleotide-overlap: number of bases overlapping [default]
#    segment-overlap: number of intervals intervals in the segments of interest overlapping annotations. A single base-pair overlap is sufficient.
#    segment-mid-overlap: number of intervals in the segments of interest overlapping at their midpoint annotations.
#    annotations-overlap: number of intervals in the annotations overlapping segments of interest. A single base-pair overlap is sufficient.
#    segment-mid-overlap: number of intervals in the annotations overlapping at their midpoint segments of interest

#data
PISOCHORE="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/d_GAT/hg19.fasta.isochore.bed.gz";
#PISOCHORE="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/d_GAT/b37.fasta.isochore.bed.gz";
#PWSP="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/d_GAT/contigs_ungapped.bed.gz";

#use as a background all the background of 1000g  used to carry out the differential analyses
#reduced by removing variants in coding regions
#PWSP="/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/ALLELESEQ_20samples_merged.vcf_sorted_clean.GAT_hg19_gencodeV16_non_cds.gz";

#VERY stringent background to test over-representation of VDR-BV in motifs
#actually I use this for everything in FIGURE 5 now: motifs, chromHMM
#basically these are all test variants WITH THE POTENTIAL of being VDR-BVs, of which only some are
#all 1kg varians under cpo3 + all VDR asb tested asym variants NOT under cpo3 (so are our recurrent strong VDR-BV enriched in VDR motifs compared to all variants under VDR?)
#PWSP="/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/GAT_BACKGROUND_ALLELESEQ_20samples_merged_INTERSECT_VDR_CPO3_peaks.bed.gz";

#The above is wrong. The following contains ~100,000 variants which were tested by AlleleSeq and are EITHER SYM or ASYM. So they have the potential to be VDR-BV
#The only problem is that many of these are not under Cpo3 peaks. So you cannot state anymore that chromHMM annotation is "OVER AND ABOVE" what we saw before
#PWSP="/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/GAT_BACKGROUND_SIMASYM_NEW_hg19.bed.gz";


if [ ! $# == 5 ]; then
        echo "Usage: `basename $0` <SEGMENT_PATH> <ANNOTATION> <BACKGROUND> <NUM_SAMPLES> <DESC>"
        echo "<SEGMENT_PATH> path to directory containing the bed file(s)"
	echo "<ANNOTATION> full path to annotation bed file"
	echo "<BACKGROUND> full path to background bed file"
	echo "<NUM_SAMPLES> number of samples to use for the randomisation (eg 50000)"
	echo "<DESC> short 1word description of the annotation used (eg wgEncodeRegDnaseClustered) - used for naming the out dir"
        exit
fi

PDATA=$1;
PANNOTATION=$2;
PWSP=$3
PRAND=$4;
DESC=$5;

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
	#              --isochore-file=${PISOCHORE}   \\
	echo "${PCODE}/gat-run.py     \\
              --ignore-segment-tracks        \\
              --segments=${FILE}             \\
              --counter=segment-overlap      \\
              --annotations=${PANNOTATION}   \\
              --workspace=${PWSP}            \\
              --isochore-file=${PISOCHORE}   \\
              --output-counts-pattern=${POUT}/${ID}.%s.oc.tsv.gz \\
              --num-samples=${PRAND}         \\
              --log=${POUT}/GAT_${DESC}_vs_${ID}.log" >> ${PDATA}/${SCRIPT};
        nice -5 qsub  -v "BASH_ENV=~/.bashrc" -e ${POUT}/GAT_${DESC}_vs_${ID}.err -o ${POUT}/GAT_${DESC}_vs_${ID}.tsv -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};  
done
