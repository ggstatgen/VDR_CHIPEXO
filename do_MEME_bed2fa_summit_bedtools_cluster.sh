#!/bin/bash

#4/10/2013

#get fasta sequences from macs2 peak summits
#this is used as an input to MEME or other motif finders
#I use the peak summit position from the best peaks in peak summit
#I add an amount of sequence left and right
#I pass it to bedtools to get the sequence

#input: summit.bed
#output: fasta intervals

#-g 2540757438 --keep-dup=all -s 40 --nomodel --extsize 26
#
#done for GATK
#
#1 - peakset beds have been concatenated
#2 - sort -k1,1V -k2,2g
#3 - bedtools slop -b 1000
#4 - bedtools merge
#
#Final "consensus" interval contains ~25000 regions and is in d_beds
#Using this as the interval file for GATK
#
#
#/net/isi-scratch/giuseppe/tools/bedtools-2.17.0/bin/bedtools slop -i all.bed  -g /net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes -b 1000 > all_1000.bed
#sort -k1,1V -k2,2g all_1000.bed > all_1000_sorted.bed
#/net/isi-scratch/giuseppe/tools/bedtools-2.17.0/bin/bedtools merge -i all_1000_sorted.bed > GM_consensus_1000.bed


if [ ! $# == 3 ]; then
        echo "Usage: `basename $0` <PATH> <LIM> <SLOP>"
        echo "<PATH> location for the summit.bed MACS2 file(s)"
	echo "<LIM> number of peaks to retain"
	echo "<SLOP> extend left and right of the summit by this amount in bp (e.g. 100)"
        exit
fi

PDATA=$1;
NPEAKS=$2;
SLOPWIN=$3;
PSCRIPTS="/net/isi-backup/giuseppe/scripts";
PBEDTOOLS="/net/isi-scratch/giuseppe/tools/bedtools-2.17.0/bin";
PCHROMSIZES="/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes";
PFASTA_M="/net/isi-scratch/giuseppe/indexes/Hsap/hg19_masked/hg19_masked.fa";
PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";

for FILE in ${PDATA}/*.bed;
        do ID=`basename ${FILE} ".bed"`;

        SCRIPT=MEME_bed2fa_summit_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};

	#extend
	echo "${PBEDTOOLS}/bedtools slop -i ${FILE}  -g ${PCHROMSIZES} -b ${SLOPWIN} > ${PDATA}/${ID}.slop${SLOPWIN}.bed" >> ${PDATA}/${SCRIPT};
	#sort
	#echo "sort -k1,1V -k2,2g ${PDATA}/${ID}.slop${SLOPWIN}.bed > ${PDATA}/${ID}.slop${SLOPWIN}.sorted.bed" >> ${PDATA}/${SCRIPT};
	#merge - maybe not needed - there will be repeated sequences, but maybe this is good?
	#echo "${PBEDTOOLS}/bedtools merge -i  ${PDATA}/${ID}.slop${SLOPWIN}.sorted.bed  > ${PDATA}/${ID}.slop${SLOPWIN}.bed" >> ${PDATA}/${SCRIPT};
	
	#ordering by signal score
	echo "cat ${PDATA}/${ID}.slop${SLOPWIN}.bed | sort -k5nr | head -n ${NPEAKS} > ${PDATA}/${ID}.best${NPEAKS}.pvalsort_summits.bed" >> ${PDATA}/${SCRIPT};
	#use bedtools to extract the fasta
	echo "${PBEDTOOLS}/bedtools getfasta -fi ${PFASTA_M} -bed  ${PDATA}/${ID}.best${NPEAKS}.pvalsort_summits.bed -fo ${PDATA}/${ID}.best${NPEAKS}.fa" >> ${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/MEME_bed2fa_summit_${ID}.err -o ${PDATA}/MEME_bed2fa_summit_${ID}.out -q newnodes.q ${PDATA}/${SCRIPT};
	#rm ${PDATA}/${ID}.slop${SLOPWIN}.bed;
	#rm ${PDATA}/${ID}.best${NPEAKS}.pvalsort_summits.bed;
        rm ${PDATA}/${SCRIPT};
done
