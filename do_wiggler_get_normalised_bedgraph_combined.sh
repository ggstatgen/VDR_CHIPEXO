#!/bin/bash

#30/5/2014
#script to run WIGGLER (also known as align2rawsignal)
#http://code.google.com/p/align2rawsignal/
#on all bam files in a directory.
#This is used to produce tracks you will then upload to UCSC browser to get normalised coverage of peaks for building images on the effect of VDT-QTL

#Ideally, you have looked at the correlation between some wigs as specified
#Then, choose the files which correlate well and run this script on them

#The program needs mappability files and so far it can ONLY DEAL with bams having mapped reads of minimum length = 20
#THEREFORE BEFORE RUNNING THIS PLEASE FILTER YOUR BAMS TO EXCLUDE MAPPED READS OF LENGTH < 20 with do_bam_filter_by_readlength_cluster.sh

#THIS WILL ONLY WORK IF THE ENVIRONMENT IN THE BASHRC IS SET SO SET IT

#OPTIONS
#----------------
#INPUT OPTIONS:
#----------------
#
#-i=<alignFname> (MANDATORY, MULTIPLE ALLOWED)
#One or more tagAlign/BAM files (replicates) as input.
#The tagAlign files can have extensions (gz,tagAlign). BAM files MUST have extension (bam,bam.gz)
#
#-s=<seqDir> (MANDATORY)
#Full path to directory containing chromosome fasta files (eg. chr1.fa ..)
#The file names MUST match the chromosome names used in the tagAlign files
#e.g: /seq/hg19
#
#-u=<uMapDir> (MANDATORY)
#Full path to directory containing binary mappability tracks.
#The directory name must be of the form [PATH]/globalmap_k<min>tok<max>
#e.g: /umap/hg19/globalmap_k20tok54
#
#----------------
#OUTPUT OPTIONS:
#----------------
#
#-o=<oFname> (OPTIONAL)
#Full path and name of output signal file.
#Set to stdout if you want to print to stdout which is also the default behavior
#Default: stdout
#
#-of=<outputFormat> (OPTIONAL)
#Output signal file format
#wiggle (wig) or bedGraph (bg) or matfile (mat)
#Default: mat
#
#-m=<localCumMapFile> (OPTIONAL)
#Calculate, for each position 'i' in the genome, the maximum number of uniquely mappable
#surrounding positions that contribute to the signal value at position 'i'.
#This is a function of the mappability of the surrounding positions, 
#tag extension/smoothing length and number of replicates.
#The local cumulative mappability is output in <maxTagsFile>.
#If -of=mat then, <localCumMapFile> can have the same name as <oFname>.
#In this case, the local cummap output is stored as a separate set of variables with prefix maxTags
#in the .mat file.
#Default: local cumMap is not output to a file
#
#-v=<logFile> (OPTIONAL)
#verbose mode.
#<logFile> Full path and name of file for logging.
#Set to stdout/stderr if you want to output logging info to stdout/stderr 
#Default: off
#
#-n=<normFLag> (OPTIONAL)
#a flag indicating whether the signal output should be normalized
#<normalization_flag> = 0,1,2,3,4,5
#Default: 5 (fold change wrt. expected value from a uniform distribution of reads)
#0: no normalization
#1: normSignal(i) = signal(i) * (1e9 / #reads)
#2: normSignal(i) = (signal(i)/winsize) * (1e9 / #reads)
#3: normSignal(i) = (signal(i)/localCumMap(i)) * (1e9 / #reads)
#4: normSignal(i) = (signal(i)/winsize) * (#total_mappable_bases / #reads)
#5: normSignal(i) = (signal(i)/localCumMap(i)) * (#total_mappable_bases / #reads)
#
#----------------
#PARAMETERS:
#----------------
#
#-l=<fragLen> (OPTIONAL, MULTIPLE ALLOWED)
#Fragment-length / 2*Tag-shift
#Default: 1 (no extension)
#Tags are shifted by floor(fragLen/2) in a 3' direction, relative to the strand of the tag
#NOTE: If a single fragLen is specified then it is applied to all tagAlign/BAM files.
#NOTE: Multiple arguments of this type are allowed. In such a case, 
#      number of fragLen arguments MUST BE == no. of Align files.
#      The ORDER of these arguments is important. 
#      e.g. The first <fragLen> is matched with the first tagAlign/BAM file and so on.
#
#-w=<smoothingWindow> (OPTIONAL)
#Smoothing window size for signal
#Default: mean(1.5*fragLen)
#
#-k=<smoothingKernel> (OPTIONAL)
#Smoothing kernel to use 
#Valid kernels (rectangular,triangular,epanechnikov,biweight,triweight,cosine,gaussian,tukey)
#Default: tukey (with taper ratio of max( 0.25 , min (0.5,max(w-mean(l),0)/(2*w)) ) if w is specified)
#
#-f=<localCumMapFilter> (OPTIONAL)
#Will nullify positions that have localCumMap <= <mappability_threshold>
#<mappability_threshold> <= 1 implies threshold is in terms of percentage localCumMap
#          e.g. 0.1 means atleast 10 percent of the positions in the smoothing window MUST be mappable 
#               > 1 implies threshold is on actual maxtags values
#          e.g. 30 means atleast 30 positions (per replicate) in the extension/smoothing window 
#               MUST be mappable 
#Default: 0.25
#
#-mm=<memory> (OPTIONAL)
#Total memory to use in GB
#Default: 2


#typical command line
#/net/isi-scratch/giuseppe/tools/align2rawsignal/bin/align2rawsignal
#-i=/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/VDR_GM06986_reads_min20.sorted.bam
#-s=/net/isi-scratch/giuseppe/indexes/Hsap/hg19/WIGGLER/seq/
#-u=/net/isi-scratch/giuseppe/indexes/Hsap/hg19/WIGGLER/globalmap_k20tok54/
#-o=/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/TEST
#-of=bg
#-v=align2rawsignal.log
#-mm=2


if [ ! $# == 2 ]; then
        echo "Usage: `basename $0` <PATH> <OUT_FORMAT>"
        echo "PATH - path for the bam files"
        echo "OUTFORMAT - one of [bg|wig|mat]"
        exit
fi

PDATA=$1;
OUTF=$2;
#the following two have been downloaded from the recommended website on the WIGGLER page
#Ideally you should run these based on sex (move the chrY data somewhere else when you run it on a female)
PMAP="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/WIGGLER/globalmap_k20tok54/";
PSEQ="/net/isi-scratch/giuseppe/indexes/Hsap/hg19/WIGGLER/seq/";
PCODE="/net/isi-scratch/giuseppe/tools/align2rawsignal/bin";


#create output dir
#POUT=${PDATA}/d_MACS2_comb_${MACSDUP}_q${MACSQ};
#mkdir ${POUT};

PREFIX=' -i=';
for FILE in ${PDATA}/*.bam;
        do
        INPUT+=${PREFIX}${FILE};
done

SCRIPT=align2raw_combo.sh;
echo '#!/bin/bash' >>${PDATA}/${SCRIPT};
echo 'source /home/giuseppe/.bashrc' >>${PDATA}/${SCRIPT};
echo "${PCODE}/align2rawsignal ${INPUT} -s=${PSEQ} -u=${PMAP} -o=${PDATA}/align2rawsignal.bg -of=${OUTF} -v=${PDATA}/align2rawsignal.log -mm=15" >>${PDATA}/${SCRIPT};

nice -5 qsub -e ${PDATA}/align2rawsignal_combo.err -o ${PDATA}/align2rawsignal_combo.out -q fgu217.q ${PDATA}/${SCRIPT};
rm ${POUT}/${SCRIPT};
