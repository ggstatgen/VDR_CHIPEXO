#!/bin/bash

#4/10/2013

#wrapper for do_MEME_output_genomic_regions_from_calls.pl
#which I hacked from the script provided by the QUEST peak caller author

#Script to create .fa files to input in MEME-chip.
#this should  create one shell script per narrowpeak file which does the following:
#sort the narrowpeak by p-value/q-value (need to decide on that)
#MACS2 forum suggestion: "remove regions with fold enrichment <5 (ordered by q-value or p-value)"
#get top X peaks (user defined)
#pass x top peaks to do_MEME_output_genomic_regions_from_calls.pl

#narropeak fields of interest
#7 signalValue - Measurement of overall (usually, average) enrichment for the region.
#8 pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
#9 qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.

#output is one fasta file per narrowPeak file. These are inputs for the MEME-suite

#I don't check if LIM is > than peak number for all files. In that case the awk script should fail and I'll get an error in the STDERR

if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <LIM> <WIN>"
        echo "<PATH> location for the .narrowPeak file(s)"
	echo "<LIM> number of peaks to retain"
	echo "<WIN> window size for peak extension (in bp, eg 100)"
        exit
fi

PDATA=$1;
NPEAKS=$2;
WIN=$3;
PSCRIPTS="/net/isi-backup/giuseppe/scripts";
PBEDTOOLS="/net/isi-scratch/giuseppe/tools/bedtools-2.17.0/bin";
PFASTA_M="/net/isi-scratch/giuseppe/indexes/Hsap/hg19_masked/hg19_masked.fa";
PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";

for FILE in ${PDATA}/*.narrowPeak;
        do ID=`basename ${FILE} ".narrowPeak"`;

        SCRIPT=MEME_bed2fa_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#ordering by p value
	#| head -n 10 | sort -k8nr
	echo "cat ${FILE} | sort -k8nr | head -n ${NPEAKS} > ${PDATA}/${ID}.best${NPEAKS}.pvalsort.narrowPeak" >> ${PDATA}/${SCRIPT};
	echo "${PSCRIPTS}/MEME_output_genomic_regions_from_calls.pl -i ${PDATA}/${ID}.best${NPEAKS}.pvalsort.narrowPeak -o ${PDATA}/${ID} -r ${PFASTA_M} -w ${WIN}" >> ${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/MEME_bed2fa_${ID}.err -o ${PDATA}/MEME_bed2fa_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done

