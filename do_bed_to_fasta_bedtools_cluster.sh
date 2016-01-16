#!/bin/bash

#4/10/2013
#script to extract the fasta sequence under peaks
#this is an alternative to do_MEME_output_genomic_regions_from_calls.pl
#that one extends the peaks and doesn't use bedtools.
#this one uses bedtools fastafrombed

#Tool:    bedtools getfasta (aka fastaFromBed)
#Version: v2.17.0
#Summary: Extract DNA sequences into a fasta file based on feature coordinates.
#
#Usage:   bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf> -fo <fasta> 
#
#Options: 
#	-fi	Input FASTA file
#	-bed	BED/GFF/VCF file of ranges to extract from -fi
#	-fo	Output file (can be FASTA or TAB-delimited)
#	-name	Use the name field for the FASTA header
#	-split	given BED12 fmt., extract and concatenate the sequencesfrom the BED "blocks" (e.g., exons)
#	-tab	Write output in TAB delimited format.
#		- Default is FASTA format.
#
#	-s	Force strandedness. If the feature occupies the antisense,
#		strand, the sequence will be reverse complemented.
#		- By default, strand information is ignored.



if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <EXT>"
        echo "<PATH> location for the interval file(s)"
	echo "<EXT> extension for the file interval (eg bed,narrowPeak)"
	echo "NOTE: change the reference before using!!"
        exit
fi

PDATA=$1;
PEXT=$2;
#location of macs2 executable, change if needed
PBEDTOOLS="/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin";
PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19_masked/hg19_masked.fa";
#PFASTA="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";
#PFASTA="/net/isi-scratch/giuseppe/indexes/Mmus/mm10/mm10.fa";
#PFASTA="/net/isi-scratch/giuseppe/indexes/Mmus/mm10_masked/mm10_masked.fa";

for FILE in ${PDATA}/*.${PEXT};
        do ID=`basename ${FILE} ".${EXT}"`;

        SCRIPT=bed2fa_${ID}.sh;

        echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo 'source activate'  >>${PDATA}/${SCRIPT};
	echo "${PBEDTOOLS}/bedtools getfasta -fi ${PFASTA} -bed ${FILE} -fo ${PDATA}/${ID}.fa" >> ${PDATA}/${SCRIPT};

        nice -5 qsub -e ${PDATA}/bed2fa_${ID}.err -q medium_jobs.q ${PDATA}/${SCRIPT};
        rm ${PDATA}/${SCRIPT};
done

