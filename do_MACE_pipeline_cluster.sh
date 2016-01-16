#!/bin/bash

#Script to run the MACE chip-exo peak caller on the VDR data
#http://dldcc-web.brc.bcm.edu/lilab/MACE/docs/html/
#input: BAM files
#output: BED files

#STEP 1
#---------Preproprocessing -  includes sequencing depth normalization, nucleotide composition bias correction, signal consolidation and noise reduction
#Replicate BAM files are separated by ','. 
#Reads mapped to forward and reverse strand will define two boundaries of binding region, so finally two wiggle files representing the coverage signal by 5' ends of reads will be produced. 
#wiggle files will be converted into bigwig format automatically if 'WigToBigWig' executable can be found in your system $PATH.:
#eg
#preprocessor.py -i CTCF_replicate1.sorted.bam,CTCF_replicate2.sorted.bam,CTCF_replicate3.sorted.bam -r hg19.chrom.sizes -o CTCF_MACE
#OPTIONS
# -i INPUT_FILE, --inputFile=INPUT_FILE
#                        Input file in BAM format. BAM file must be sorted and
#                        indexed using samTools. Replicates separated by
#                        comma(',') e.g. "-i rep1.bam,rep2.bam,rep3.bam"
#  -r CHROMSIZE, --chromSize=CHROMSIZE
#                        Chromosome size file. Tab or space separated text file
#                        with 2 columns: first column is chromosome name,
#                        second column is size of the chromosome.
#  -o OUTPUT_PREFIX, --outPrefix=OUTPUT_PREFIX
#                        Prefix of output wig files(s). "Prefix_Forward.wig"
#                        and "Prefix_Reverse.wig" will be generated
#  -w WORD_SIZE, --kmerSize=WORD_SIZE
#                        Kmer size [6,10] to correct nucleotide composition
#                        bias. Kmer size will significantly affects running
#                        speed. For each read in BAM file, MACE needs to lookup
#                        table that has 4^Kmer items. Read_length > word_size*2
#                        default=6
#  -b BIN, --bin=BIN     Chromosome chunk size. Each chomosome will be cut into
#                        small chunks of this size. Decrease chunk size will
#                        save more RAM. default=100000 (bp)
#  -d REFREADN, --depth=REFREADN
#                        Reference reads count (default = 10 million).
#                        Sequencing depth will be normailzed to this number, so
#                        that wig files are comparable between replicates.
#  -q QUAL_CUT, --qCut=QUAL_CUT
#                        phred scaled mapping quality threshhold to determine
#                        "uniqueness" of alignments. default=30

PP_OPTS="-w 6 -d 20000000";

#STEP 2
#-----------Border detection and border pairing
#eg mace.py -s hg19.chrom.sizes -f CTCF_MACE_Forward.bw -r CTCF_MACE_Reverse.bw -o CTCF_MACE
#OPTIONS
#  -f FORWARD_BW, --forward=FORWARD_BW
#                        BigWig format file containing coverage calcualted from
#                        reads mapped to *forward* strand.
#  -r REVERSE_BW, --reverse=REVERSE_BW
#                        BigWig format file containing coverage calcualted from
#                        reads mapped to *reverse* strand.
#  -s CHROMSIZE, --chromSize=CHROMSIZE
#                        Chromosome size file. Tab or space separated text file
#                        with 2 columns: first column contains chromosome name,
#                        second column contains chromosome size. Example:chr1
#                        249250621 <NewLine> chr2        243199373 <NewLine>
#                        chr3        198022430 <NewLine> ...
#  -o OUTPUT_PREFIX, --out-prefix=OUTPUT_PREFIX
#                        Prefix of output files. NOTE: if 'prefix.border.bed'
#                        exists and was non-empty, peak calling step will be
#                        skipped! So if you want to rerun mace.py from scratch,
#                        use different 'prefix' or delete old
#                        'prefix.border.bed' before starting.
#  -p PVALUE_CUTOFF, --pvalue=PVALUE_CUTOFF
#                        Pvalue cutoff for border detection and subsequent
#                        border pairing. default=0.05
#  -m MAX_DISTANCE, --max-dist=MAX_DISTANCE
#                        Maximum distance allowed for border pairing.
#                        default=100
#  -e FUZZY_SIZE, --fz-window=FUZZY_SIZE
#                        Peaks located closely within this window will be
#                        merged. default=5 (bp)
#  -w WINDOW_SIZE, --bg-window=WINDOW_SIZE
#                        Background window size used to determine background
#                        signal level. default=100 (bp)
#  -n SIGNAL_FOLD, --fold=SIGNAL_FOLD
#                        Minmum coverage signal used to build model (i.e.
#                        estimate optimal peak pair size). default=2.0
M_OPTS="-n 1.2";


if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <CHROM_INFO>"
        echo "<PATH> data path for the bam files"
        echo "CHROM_INFO - file of chromosome sizes"
        exit
fi

PDATA=$1;
PCHROM=$2
#location of MACE executables, change if needed
PCODE="/home/giuseppe/src/pypa-virtualenv-3a798b3/ve/bin";


n=$RANDOM;
POUT=${PDATA}/d_MACE_$n;
mkdir ${POUT};
echo "pre: ${PP_OPTS}  post: ${M_OPTS}" >> ${POUT}/setup.log;

for FILE in ${PDATA}/*.bam;
	do ID=`basename ${FILE} ".bam"`;

        SCRIPT=MACE_1_${ID}.sh;

        echo '#!/bin/bash' >>${POUT}/${SCRIPT};
        echo '' >>${POUT}/${SCRIPT};
	echo 'source activate' >> ${POUT}/${SCRIPT};
	echo "${PCODE}/preprocessor.py -i ${FILE} -r ${PCHROM} -o ${POUT}/MACE_${ID} ${PP_OPTS}" >>${POUT}/${SCRIPT};
	#wig to bigwig?	
	#mace_GM07029_Forward.wig
	#mace_GM07029_Reverse.wig
	#mace.py -s hg19.chrom.sizes -f CTCF_MACE_Forward.bw -r CTCF_MACE_Reverse.bw -o CTCF_MACE
	#potential opts: -p (pvalue) 
	echo "${PCODE}/mace.py -s ${PCHROM} -f ${POUT}/MACE_${ID}_Forward.bw -r ${POUT}/MACE_${ID}_Reverse.bw  -o ${POUT}/MACE_${ID} ${M_OPTS}" >>${POUT}/${SCRIPT};
	#output:  prefix.border_pair.bed

        nice -5 qsub -e ${POUT}/MACE_${ID}.err -o ${POUT}/MACE_${ID}.out -q newnodes.q ${POUT}/${SCRIPT};
        rm ${POUT}/${SCRIPT}; 
done

