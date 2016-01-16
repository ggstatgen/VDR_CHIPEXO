#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

#11/3/2013
#testing the idr stuff
#this is used to do the idr thresholding on a set of pseudo replicates of a peak file
#it follows this
#https://sites.google.com/site/anshulkundaje/projects/idr#TOC-IDR-ANALYSIS-ON-SELF-PSEUDOREPLICATES
#this uses other scripts:
#1 do_bam_to_tagalign.sh
#2 do_random_split_peakfile.sh
#3 idr_count_peaks.sh

my $input_id;
my $datadir;
if ($#ARGV < 0) {
     print "USAGE: _idr_pseudorep.pl <INFILE>\n";
     print "<INFILE> .bam file\n";
     exit 1;
}
my $infile=$ARGV[0];

if($infile =~ /(VDR\_)(\w{2}\d{5})(\_Peconic\d{5})(\_trimmed)(.bam)$/){  #eg VDR_GM19213_Peconic20324_trimmed.bam
	$input_id = $2;
}else{
	print "Input file name not recognised. Modify the regular expression in the script.\n";
	exit -1;
}

my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
if($directory eq "\.\/"){ 
	$datadir = $input_id . '_idr_thresholding';
}else{
	$datadir =  $directory . "\/" . $input_id . '_idr_thresholding';
}
system "mkdir $datadir";
#---------
#temp data
#---------
my $infile_tagAlign    = $datadir . "\/" . $basename . '.tagAlign';
#my $infile_tagAlign_gz = $datadir . "\/" . $basename . '.tagAlign.gz';
my $infile_tagAlign_gz =  $basename . '.tagAlign.gz';
my $infile_tagAlign_pr1 = $basename . '.pr1.tagAlign.gz';
my $infile_tagAlign_pr2 = $basename . '.pr2.tagAlign.gz';
my $macs2_basename         = $datadir . "\/" . $basename  . '_macs2';  #macs2 outputs in these
my $macs2_pr1_basename      = $datadir . "\/" . $basename  . '_macs2_pr1';  #macs2 outputs here for pseudoreplicate 1
my $macs2_pr2_basename      = $datadir . "\/" . $basename  . '_macs2_pr2';  #macs2 outputs here for pseudoreplicate 1
my $macs2_encodePeak_file     = $macs2_basename     . '_peaks.encodePeak';;
my $macs2_pr1_encodePeak_file = $macs2_pr1_basename . '_peaks.encodePeak';
my $macs2_pr2_encodePeak_file = $macs2_pr2_basename . '_peaks.encodePeak';
my $macs2_regionPeak_file        = $macs2_basename     . '.regionPeak';
my $macs2_pr1_regionPeak_file    = $macs2_pr1_basename . '.regionPeak';
my $macs2_pr2_regionPeak_file    = $macs2_pr2_basename . '.regionPeak';
my $macs2_regionPeak_gz_file     = $macs2_basename     . '.regionPeak.gz';
my $macs2_pr1_regionPeak_gz_file = $macs2_pr1_basename . '.regionPeak.gz';
my $macs2_pr2_regionPeak_gz_file = $macs2_pr2_basename . '.regionPeak.gz';
my $idr_basename = $datadir . "\/" . 'pr1VSpr2';
my $idr_plot_basename = $datadir . "\/"  . 'pr1VSpr2_plot';
my $idr_overlap_file = $idr_basename . '-overlapped-peaks.txt';


my $SCRIPT_PATH = '/net/isi-scratch/giuseppe/scripts';
my $TOOL_PATH = '/net/isi-scratch/giuseppe/tools';
#my $CODE_PATH = '/home/giuseppe/local/bin';
my $R_CODE = '/net/isi-cgat/ifs/apps/apps/R-2.14.1/bin';
#program paths, change as needed
my $PICARD = $TOOL_PATH . '/picard-tools-1.84/picard-tools-1.84';
my $SAMTOOLS = $TOOL_PATH . '/samtools';
my $BEDTOOLS = $TOOL_PATH . '/bedtools-2.17.0/bin';
my $SPP = $TOOL_PATH . '/phantompeakqualtools';
my $IDR = $TOOL_PATH . '/idrCode';
my $MAPPABILITY = 2540757438;


#bam -> tagalign
print "-----------------------\n";
print "INFO - create gzipped tagAlign from initial Bam file\n";
print "-----------------------\n";
system "$SCRIPT_PATH/do_bam_to_tagalign.sh $infile";

#Call peaks  (REMEMBER TO USE RELAXED THRESHOLDS AND TRY TO CALL 150k to 300k peaks even if most of them are noise)
#system "$R_CODE/Rscript $SPP/run_spp.R -c=$infile_tagAlign_gz -i=controlSampleRep0.tagAlign.gz -npeak=300000 -odir=/peaks/reps -savr -savp -rf -out=/stats/phantomPeakStatsReps.tab";
#or alternatively MACS2 called as follows
#system "macs2 callpeak -t $infile -n $macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 -p 1e-3 --to-large";
#or if you want to use the tagAlign, do
print "-----------------------\n";
print "INFO - MACS2: Call peaks full file\n";
print "-----------------------\n";
system "macs2 callpeak -t $infile_tagAlign_gz -f BED -n $macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 -p 1e-3 --to-large";

#MACS2 will output peaks in various formats (bed, xls and encodePeak). You want to use the *encodePeak files. 
#These are in the standard UCSC/ENCODE narrowPeak format. Sort the *encodPeak files from best to worst using the -log10(pvalue) column i.e. column 8. 
#Then truncate the number of peaks to the top 100k-125k. Using more than this simply increases the running time of the IDR analysis with no advantage.
#Infact using more peaks with MACS2 can cause problems with the IDR model because MACS2 seems to produce strange highly correlated peak scores for very weak and noisy detections.
#This can confuse the IDR model.

#sort MACS2 output encodePeak file
print "-----------------------\n";
print "INFO - sort and gzip MACS2 output regionPeak file\n";
print "-----------------------\n";
system "sort -k 8nr,8nr $macs2_encodePeak_file | head -n 100000 | gzip -c > $macs2_regionPeak_file";

#Randomly split the mapped reads in the tagAlign file into 2 parts (pseudoReplicates) with approximately equal number of reads
#use script
print "-----------------------\n";
print "INFO - create two random subsets of tagAlign file\n";
print "-----------------------\n";
system "$SCRIPT_PATH/do_random_split_peakfile.sh $infile_tagAlign_gz";

#you should now have two tagAlign_gz file
#use macs2 to call peaks as done above
#each one becomes a thread obviously
print "-----------------------\n";
print "INFO - Call peaks pseudoreplicate 1\n";
print "-----------------------\n";
system "macs2 callpeak -t $infile_tagAlign_pr1 -f BED -n $macs2_pr1_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 -p 1e-3 --to-large";
print "-----------------------\n";
print "INFO - Call peaks pseudoreplicate 2\n";
print "-----------------------\n";
system "macs2 callpeak -t $infile_tagAlign_pr2 -f BED -n $macs2_pr2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 -p 1e-3 --to-large";


print "-----------------------\n";
print "INFO - sort and gzip MACS2 output regionPeak file (pseudorep 1)\n";
print "-----------------------\n";
#system "sort -k 8nr,8nr $macs2_pr1_encodePeak_file | head -n 100000 | gzip -c > $macs2_pr1_regionPeak_file";
system "sort -k 8nr,8nr $macs2_pr1_encodePeak_file | head -n 100000 > $macs2_pr1_regionPeak_file";

print "-----------------------\n";
print "INFO - sort and gzip MACS2 output regionPeak file (pseudorep 2 )\n";
print "-----------------------\n";
#system "sort -k 8nr,8nr $macs2_pr2_encodePeak_file | head -n 100000 | gzip -c > $macs2_pr2_regionPeak_file";
system "sort -k 8nr,8nr $macs2_pr2_encodePeak_file | head -n 100000 > $macs2_pr2_regionPeak_file";

###
#IDR Analysis
###
#For IDR analysis, we always compare a pair of peak files
#At this point you should have
#1 VDR_GM19235_Peconic20333_trimmed_macs2.regionPeak.gz
#This is useless because you have no other file to compare it to (or maybe against pr1 and then pr2?)
#2 the two subsamples (PSEUDO REPLICATES)
#VDR_GM19235_Peconic20333_trimmed_macs2_pr1.regionPeak.gz
#VDR_GM19235_Peconic20333_trimmed_macs2_pr2.regionPeak.gz


#MACS2: If you are using peaks based on the MACS2 peak caller, then use p.value as the ranking measure.
print "-----------------------\n";
print "INFO - IDR run on pseudoreplicates\n";
print "-----------------------\n";
#usage
#Rscript batch-consistency-analysis.r [peakfile1] [peakfile2] [peak.half.width] [outfile.prefix] [min.overlap.ratio] [is.broadpeak] [ranking.measure]
#Rscript batch-consistency-analysis.r [peakfile1] [peakfile2] -1 [outfile.prefix] 0 F p.value
system "$R_CODE/Rscript batch-consistency-analysis.r $macs2_pr1_regionPeak_file $macs2_pr2_regionPeak_file -1 $idr_basename 0 F p.value";

print "-----------------------\n";
print "INFO - IDR consistency plot\n";
print "-----------------------\n";
#Rscript batch-consistency-plot.r [npairs] [output.prefix] [input.file.prefix1] [input.file.prefix2] [input.file.prefix3] ....
#[n.pairs] is the number of pairs of replicates that you want to plot on the same plot
#e.g. 1 or 3 or ...
#[output.prefix] is a prefix that will be used to name output data from this analysis. 
#NOT TO BE CONFUSED with [outfile.prefix] in batch-consistency-analysis.r
#The prefix must also include the PATH to the directory where you want to store the output data.
#e.g. /consistency/plots/chipSampleAllReps
system "$R_CODE/Rscript batch-consistency-plot.r 1 $idr_plot_basename  $idr_basename";
#===================================================
#GETTING NUMBER OF PEAKS THAT PASS AN IDR THRESHOLD
#===================================================
#For each pairwise analysis, we have a *overlapped-peaks.txt file
#
#The last column (Column 11) of the overlapped-peaks.txt file has the global IDR score for each pair of overlapping peaks
#To get the number of peaks that pass an IDR threshold of T (e.g. 0.01) you simply find the number of lines that have a global IDR score <= T
#For self-consistency analysis of datasets with shallow sequencing depths, you can use an IDR threshold as relaxed as 0.1 if you start with < 100K pre-IDR peaks.
my $peak_number = `./idr_count_peaks.sh  $idr_overlap_file`;
print $peak_number, "\n";