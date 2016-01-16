#!/bin/bash
#uses cutadapt to remove adapter sequences from reads.
# this is using adapter sequencing hardcoded here and found to be enriched in the VDR dataset

#28/3
#I reduced the number of adapters to those that were overrepresented
#I also added reverse complement sequences

#usage cutadapt -e ERROR-RATE -a ADAPTER-SEQUENCE input.fastq > output.fastq
#e parameter: attention, now you're allowing sequencing error to be detected in the adapter sequence. Correct?

#06/12/12 added other sequences taken from the fastqc list
# these were enriched in the document I sent to Ram etc, together with the truseq adapters
# Illumina Multiplexing PCR Primer 2.01			GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
# Illumina PCR Primer Index 1						CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC
# Illumina PCR Primer Index 2						CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC
# Illumina PCR Primer Index 3						CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC
# Illumina PCR Primer Index 4						CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC
# Illumina PCR Primer Index 5						CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC
# Illumina PCR Primer Index 6						CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC
# Illumina PCR Primer Index 7						CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC
# Illumina PCR Primer Index 8						CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC
# Illumina PCR Primer Index 9						CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC
# Illumina PCR Primer Index 10					CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC
# Illumina PCR Primer Index 11					CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC
# Illumina PCR Primer Index 12					CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC
# added them at the beginning of the list

#cutadapt options
#    -b ADAPTER, --anywhere=ADAPTER
#                        Sequence of an adapter that was ligated to the 5' or
#                        3' end. If the adapter is found within the read or
#                        overlapping the 3' end of the read, the behavior is
#                        the same as for the -a option. If the adapter overlaps
#                        the 5' end (beginning of the read), the initial
#                        portion of the read matching the adapter is trimmed,
#                        but anything that follows is kept.

#    -e ERROR_RATE, --error-rate=ERROR_RATE
#                        Maximum allowed error rate (no. of errors divided by
#                        the length of the matching region) (default: 0.1)

#    -O LENGTH, --overlap=LENGTH
#                        Minimum overlap length. If the overlap between the
#                        read and the adapter is shorter than LENGTH, the read
#                        is not modified.This reduces the no. of bases trimmed
#                        purely due to short random adapter matches (default:
#                        3).

#    -m LENGTH, --minimum-length=LENGTH
#                        Discard trimmed reads that are shorter than LENGTH.
#                        Reads that are too short even before adapter removal
#                        are also discarded. In colorspace, an initial primer
#                        is not counted (default: 0).

#    -q CUTOFF, --quality-cutoff=CUTOFF
#                        Trim low-quality ends from reads before adapter
#                        removal. The algorithm is the same as the one used by
#                        BWA (Subtract CUTOFF from all qualities; compute
#                        partial sums from all indices to the end of the
#                        sequence; cut sequence at the index at which the sum
#                        is minimal) (default: 0)


if [ ! $# == 1 ]; then
        echo "Usage: `basename $0` <PATH_DATA>"
        echo "PATH_DATA - directory containing the input fastq.gz"
        exit
fi

PDATA=$1;
PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";
#my new version of cutadapt is there

for FILE in ${PDATA}/*.fastq.gz;
        do ID=`echo ${FILE} | egrep -o "GM[0-9]*"`;

	BASEFILE=`basename ${FILE} ".fastq.gz"`;
	
	SCRIPT=${ID}.sh;
	echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
        echo '' >>${PDATA}/${SCRIPT};
	echo "${PCODE}/cutadapt -b GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG  -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -b CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -b CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -b CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -b CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -b CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -b CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -b CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -b CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -b CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -b CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -b CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -O 9 -m 20 -e 0.15 ${FILE} -o ${PDATA}/${BASEFILE}_trimmed.fastq.gz" >>${PDATA}/${SCRIPT};
        nice -5 qsub -e ${PDATA}/${ID}.err -o ${PDATA}/${ID}.out -q medium_jobs.q ${PDATA}/${SCRIPT};
	rm ${PDATA}/${SCRIPT}; 
done
