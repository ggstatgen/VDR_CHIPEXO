#!/bin/bash

#script to run peakzilla on all the VDR samples, using them as replicates
#Usage: python peakzilla.py [OPTIONS] chip.bed control.bed > results.tsv
#Options:
#  -h, --help            show this help message and exit
#  -m N_MODEL_PEAKS, --model_peaks=N_MODEL_PEAKS
#                        number of most highly enriched regions used to
#                        estimate peak size: default = 200
#  -c ENRICHMENT_CUTOFF, --enrichment_cutoff=ENRICHMENT_CUTOFF
#                        minimum cutoff for fold enrichment: default = 2
#  -s SCORE_CUTOFF, --score_cutoff=SCORE_CUTOFF
#                        minimum cutoff for peak score: default = 1
#  -f FRAGMENT_SIZE, --fragment_size=FRAGMENT_SIZE
#                        manually set fragment size in bp: default = estimate
#                        from data
#  -e, --gaussian        use empirical model estimate instead of gaussian
#  -p, --bedpe           input is paired end and in BEDPE format
#  -l LOG, --log=LOG     directory/filename to store log file to: default =
#                        log.txt
#  -n, --negative        write negative peaks to negative_peaks.tsv


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH>"
#	echo "Usage: `basename $0` <PATH> <N_MODEL_PEAKS> <ENRICHMENT_CUTOFF>"
        echo "PATH - path for the bed files"
        #echo "N_MODEL_PEAKS (default=200)"
        #echo "ENRICHMENT_CUTOFF (default=2)"
        exit
fi

PDATA=$1;
#MODEL_PEAKS=$2;
#ENRICH=$3;
PCODE="/net/isi-scratch/giuseppe/tools/peakzilla-master";

#create output dir
#POUT=${PDATA}/d_peakzilla_m${MODEL_PEAKS}_c${ENRICH};
#mkdir ${POUT};

for FILE in ${PDATA}/*.bed;
        do ID=`basename ${FILE} ".bed"`;

        #SCRIPT=pz_m${MODEL_PEAKS}_c${ENRICH}_${ID}.sh;
	SCRIPT=pz_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

        echo "source activate" >> ${PDATA}/${SCRIPT};
        echo "python ${PCODE}/peakzilla.py --fragment_size=30  -l ${PDATA}/peakzilla_${ID}.log ${FILE}" >> ${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/peakzilla_${ID}.err -o ${PDATA}/peakzilla_${ID}.tsv -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done

