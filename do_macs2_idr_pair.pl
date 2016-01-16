#!/usr/bin/perl
use strict;
use warnings;
use threads;
use Getopt::Long;
use File::Basename;
#18/4/2013
#This calls peaks with MACS 2 and filters them based on IDR BETWEEN two different files.

#input: 2 bam files, input 1, input 2
#output: IDR bed

my $infile1; 
my $infile2;
my $datadir;
my $bamfile1; 
my $bamfile2; 
my $input_id1; 
my $input_id2;
my $directory;
my $IDR_THRESHOLD = 0.01;
my $MAPPABILITY = 2540757438;
#my $MAPPABILITY = 1870000000; #mouse

my $SCRIPT_PATH  = '/net/isi-backup/giuseppe/scripts';
my $TOOL_PATH    = '/net/isi-scratch/giuseppe/tools';
my $R_CODE       = '/net/isi-scratch/giuseppe/tools/R-3.1.0/bin';
my $IDR          = $TOOL_PATH . '/idrCode'; 

GetOptions(
	'i1=s'  => \$infile1,
	'i2=s'  => \$infile2
);
if(!$infile1){
     print "USAGE: do_macs2_idr_pair.pl -i1=<INFILE1> -i2=<INFILE2>\n";
     print "<INFILE1> bam file 1\n";
     print "<INFILE2> bam file 2\n";
     exit 1;
}
if(!$infile2){
     print "USAGE: do_macs2_idr_pair.pl -i1=<INFILE1> -i2=<INFILE2>\n";
     print "<INFILE1> bam file 1\n";
     print "<INFILE2> bam file 2\n";
     exit 1;
}

#if($infile1 =~ /(VDR\_)(\w{2}\d{5})(\_Peconic\d{5}).*/){  #eg VDR_GM19213_Peconic20324_trimmed.fastq or VDR_GM19214_Peconic20318_trimmed.sorted.bam
#if($infile1 =~ /(VDR\_)(\w{2}\d{5})_reads.*/){ #VDR_GM12345)_reads.bam
#	$input_id1 = $2;
#}else{
#	print "Input file name 1 not recognised. Modify the regular expression in the script.\n";
#	exit -1;
#}
##if($infile2 =~ /(VDR\_)(\w{2}\d{5})(\_Peconic\d{5}).*/){
#if($infile2 =~ /(VDR\_)(\w{2}\d{5})_reads.*/){
#	$input_id2 = $2;
#}else{
#	print "Input file name 2 not recognised. Modify the regular expression in the script.\n";
#	exit -1;
#}
#GIU OCT 2013
$input_id1 = "file1";
$input_id2 = "file2";

($bamfile1, $directory) = fileparse($infile1);
($bamfile2, $directory) = fileparse($infile2);
my $basename1 = $bamfile1; 
my $basename2 = $bamfile2; 
$basename1 =~ s/(.*)\..*/$1/; $basename2 =~ s/(.*)\..*/$1/; #remove extension
#I want all the outputs in a directory whose format will be $PATH/IDR_$input_id1_$input_id2
if($directory eq "\.\/"){ 
	$datadir = 'IDR_MACS2_' . $input_id1 . '_' . $input_id2;
}else{
	$datadir =  $directory . 'IDR_MACS2_' . $input_id1 . '_' . $input_id2;
}
system "mkdir $datadir";

#file 1
my $infile_tagAlign_1_gz       = $datadir . "\/" . $basename1 . '.tagAlign.gz';
my $macs2_basename_1           = $datadir . "\/" . $basename1 . '_macs2'; 
my $macs2_encodePeak_file_1    = $macs2_basename_1 . '_peaks.narrowPeak';
my $macs2_regionPeak_file_1    = $macs2_basename_1 . '.regionPeak';
my $macs2_regionPeak_gz_file_1 = $macs2_basename_1 . '.regionPeak.gz';
#file 2
my $infile_tagAlign_2_gz       = $datadir . "\/" . $basename2 . '.tagAlign.gz';
my $macs2_basename_2           = $datadir . "\/" . $basename2  . '_macs2'; 
my $macs2_encodePeak_file_2    = $macs2_basename_2 . '_peaks.narrowPeak';
my $macs2_regionPeak_file_2    = $macs2_basename_2 . '.regionPeak';
my $macs2_regionPeak_gz_file_2 = $macs2_basename_2 . '.regionPeak.gz';

my $idr_basename                = $datadir               . "\/" . 'IDR_MACS2_' . $input_id1 . '_vs_' . $input_id2;;
my $idr_plot_basename           = $datadir               . "\/" . 'IDR_MACS2_' . $input_id1 . '_vs_' . $input_id2 . '_IDR_plot';
my $idr_overlap_file            = $idr_basename          . '-overlapped-peaks.txt';
my $idr_overlap_file_thrs       = $idr_basename          . '-overlapped-peaks-thrs.txt';
my $idr_encodepeak_file         = $idr_overlap_file_thrs . '.npk';
#final outputs
my $idr_peakset_1                     = $datadir . "\/" . $basename1 . '_idr_out.npk';
my $idr_peakset_2                     = $datadir . "\/" . $basename2 . '_idr_out.npk';
my $output_file_1                     = $datadir . "\/" . $basename1 . '_idr_out.bed';
my $output_file_2                     = $datadir . "\/" . $basename2 . '_idr_out.bed';
my $output_file_union                 = $datadir . "\/" . 'overlap_idr_out.bed';

$directory  = substr($directory, 0, -1);
#---------------
#bam -> tagalign
#---------------
#thread this
print "-----------------------\n";
print "INFO - bam -> tagalign\n";
print "-----------------------\n";
my $thread_bam_to_tagalign1  = threads->create( \&bam_to_tagalign, $infile1);
my $thread_bam_to_tagalign2  = threads->create( \&bam_to_tagalign, $infile2);
$thread_bam_to_tagalign1->join;
$thread_bam_to_tagalign2->join;
if (my $err = $thread_bam_to_tagalign1->error()) { warn("Thread bam_to_tagalign() error: $err\n"); } 
if (my $err = $thread_bam_to_tagalign2->error()) { warn("Thread bam_to_tagalign() error: $err\n"); } 
#-------------------------------
# call peaks on file pair with MACS2
#-------------------------------
print "-----------------------\n";
print "INFO - MACS2: Call peaks\n";
print "-----------------------\n";
my $thread_macs2_file1  = threads->create( \&run_macs2,  $infile_tagAlign_1_gz, $macs2_basename_1);
my $thread_macs2_file2  = threads->create( \&run_macs2,  $infile_tagAlign_2_gz, $macs2_basename_2);
$thread_macs2_file1->join;
$thread_macs2_file2->join;
if (my $err = $thread_macs2_file1->error()) { warn("Thread run_macs2() error: $err\n"); } 
if (my $err = $thread_macs2_file2->error()) { warn("Thread run_macs2() error: $err\n"); } 
#---------------------------------
#sort MACS2 output encodePeak file 1 and 2
#---------------------------------
print "-----------------------\n";
print "INFO - sort peaksets and get top 100,000 regions\n";
print "-----------------------\n";
my $thread_sort_file1  = threads->create( \&sort_peakset,  $macs2_encodePeak_file_1, $macs2_regionPeak_file_1);
my $thread_sort_file2  = threads->create( \&sort_peakset,  $macs2_encodePeak_file_2, $macs2_regionPeak_file_2);
$thread_sort_file1->join;
$thread_sort_file2->join;
if (my $err = $thread_sort_file1->error()) { warn("Thread sort_pseudorep() error: $err\n"); } 
if (my $err = $thread_sort_file2->error()) { warn("Thread sort_pseudorep() error: $err\n"); }
#------------
#IDR Analysis
#------------
print "-----------------------\n";
print "INFO - IDR run on peakset pair\n";
print "-----------------------\n";
system "$R_CODE/Rscript $IDR/batch-consistency-analysis.r $macs2_regionPeak_file_1 $macs2_regionPeak_file_2 -1 $idr_basename 0 F p.value";
system "$R_CODE/Rscript $IDR/batch-consistency-plot.r 1 $idr_plot_basename  $idr_basename";
my $IDR_peaks = `$SCRIPT_PATH/idr_get_peaks.sh  $idr_overlap_file $IDR_THRESHOLD | wc -l`; #get the number of peak which pass the threshold
print $IDR_peaks;
chomp $IDR_peaks;
#----------------
#IDR - get peak set in original macs call which pass IDR threshold
#----------------
#MACS2: Sort peaks by p.value before selecting the top N peaks from each file peak call
system "cat $macs2_encodePeak_file_1 | sort -k8nr,8nr | head -n $IDR_peaks > $idr_peakset_1";
system "cat $macs2_encodePeak_file_2 | sort -k8nr,8nr | head -n $IDR_peaks > $idr_peakset_2";
#get beds because igv does not seem to be able to read narrowPeak files
system "$SCRIPT_PATH/narrowPeak_to_bed.pl -i=$idr_peakset_1 > $output_file_1";
system "$SCRIPT_PATH/narrowPeak_to_bed.pl -i=$idr_peakset_2 > $output_file_2";
#----------------
#IDR - get candidate overlapping peak list - get subset of -overlapped-peaks.txt that pass $IDR_THRESHOLD ($full_idr_overlap_file_thrs)
#----------------------
system "$SCRIPT_PATH/idr_get_peaks.sh $idr_overlap_file $IDR_THRESHOLD > $idr_overlap_file_thrs";
system "$SCRIPT_PATH/idrOverlap2npk.sh $idr_overlap_file_thrs";
#anshul says using this will give coarser peaks than those captured by the peak caller though
#get beds because igv does not seem to be able to read narrowPeak files
system "$SCRIPT_PATH/narrowPeak_to_bed.pl -i=$idr_encodepeak_file > $output_file_union";
#clean
unlink $infile_tagAlign_1_gz; 
unlink $infile_tagAlign_2_gz;
system "rm -rf $macs2_basename_1*"; #remove macs2 output
system "rm -rf $macs2_basename_2*"; #remove macs2 output
system "rm -rf $idr_basename*"; #remove IDR output
print "FINISHED", "\n";
#-----------
#subs
#-----------
sub bam_to_tagalign{
	my ($bam_input) = @_;
	system "$SCRIPT_PATH/do_bam_to_tagalign.sh $bam_input $datadir";
	if ( $? == -1 ){
    	print "bam_to_tagalign(): problem with output: $!\n";
    	exit -1;
	}
	return;
}
#---------------------------------------------------------------------------------
sub run_macs2{
	my ( $input_file, $macs2_basename ) = @_;
	#use the following from IDR:
	#system "macs2 callpeak -t $input_file -f BED -n $macs2_basename -g $MAPPABILITY --keep-dup=1 --nomodel --extsize 26 -p 1e-3 --verbose 1";
	#system "macs2 callpeak -t $input_file -f BED -n $macs2_basename -g $MAPPABILITY --keep-dup=1 --nomodel -p 1e-3 --extsize 26 --verbose 1";
	system "macs2 callpeak -t $input_file -f BED -n $macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel -p 1e-3 --extsize 116 --verbose 1";
	if ( $? == -1 ){
		print "MACS2: problem with output: $!\n";
        exit -1;
    }
    return;	
}
#----------------------------------------------------------------------------------
sub sort_peakset{
	my ( $encodePeak_file, $regionPeak_file ) = @_;
	system "sort -k 8nr,8nr $encodePeak_file | head -n 100000 > $regionPeak_file";
	return;
}
