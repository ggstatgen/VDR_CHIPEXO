#!/usr/bin/perl
use strict;
use warnings;
use threads;
use Getopt::Long;
use File::Basename;
use Statistics::Descriptive;
use List::Util qw(sum);

#21/3/2013
#This is obtained from saturation_test_final.pl and FRiP_test.pl
#1 does  saturation_analysis, median enrichment plots, FRiP
#2 for the saturation analysis, it does genetrack(+postproc), GPS, MACS 1.4, MACS2 (+IDR)
#does the full postprocessing IDR pipeline (standalone in the _idr_pseudorep.pl file) after each MACS2 peak calling
#3 gets median enrichment
#4 gets FRiP metric

#IDR is very VERY slow and the other threads will all be waiting for the idr ones
#also Anshul suggests not to use on MACS14
#should probably also try it with GEM
#should probably parallelise the random subsets, for now they happen in sequence and all output for one iteration is deleted after output peak number is collected

#INPUTS
#1 fastq file
#2 bin size in number of mapped reads

#OUTPUTS
#TSV FILE  <samplenumber> <uniquely mapped reads> [<number of peaks> <number of peaks/total peaks for full file>] <marginal_fe> <FRiP>
#TSV file with full peak signals per iteration (rows: iterations, columns: peaks found)

my $infile;my $INITIAL_BIN;
my $input_id; my $datadir; # to customise randomised outputs
my $SCRIPT_PATH     = '/net/isi-backup/giuseppe/scripts';
my $TOOL_PATH       = '/net/isi-scratch/giuseppe/tools';
my $INDEX_PATH      = '/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19_bowtie';
my $CHROM_SIZE_PATH = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes';
my $GENOME_PATH     = '/net/isi-scratch/giuseppe/indexes/Hsap/chromosomes_hg19';
#my $R_CODE          = '/net/isi-cgat/ifs/apps/apps/R-2.14.1/bin';
my $R_CODE          = '/net/isi-scratch/giuseppe/tools/R-3.0.1/bin';
my $TMP_DIR         = '/tmp';

#program paths, change as needed
my $BOWTIE    = $TOOL_PATH . '/bowtie-0.12.9/bowtie';
my $PICARD    = $TOOL_PATH . '/picard-tools-1.102';
my $SAMTOOLS  = $TOOL_PATH . '/samtools';
my $BEDTOOLS  = $TOOL_PATH . '/bedtools-2.17.0/bin';
my $IDR       = $TOOL_PATH . '/idrCode'; #now using idr in giuseppe/scripts
my $GEM       = $TOOL_PATH . '/gem';
my $GENETRACK = $TOOL_PATH . '/chipexo-master/genetrack/genetrack.py';
#my $SPP = $TOOL_PATH . '/phantompeakqualtools';

#human- and dataset-specific
my $IDR_THRESHOLD = 0.01;
my $SLOP = 1000; #how many bp to extend the bed intervals left and right of each peak before FRiP calculation
my $MAPPABILITY = 2540757438; #for Hsap (40bp) used by Macs2 and GEM
#macs2 needs to be ran with very loose thresholds when using it with IDR
my $MACS14_OPTS    = "--keep-dup=1 --nomodel --shiftsize 13 --verbose 1";
my $MACS2_OPTS     = "--keep-dup=1 -q 0.01 --nomodel --extsize 26 --verbose 1";
my $MACS2_IDR_OPTS = "--keep-dup=1 -p 1e-3 --nomodel --extsize 26 --to-large --verbose 1";
my $GEM_OPTS       = "--genome $GENOME_PATH --k_min 3 --k_max 18 -top 0 --outBED";
my @HUMAN_CHRS = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                  'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 
                  'chr20', 'chr21', 'chr22', 'chrX', 'chrY'); 
GetOptions(
	'i=s'  => \$infile,
	'bin=f'    => \$INITIAL_BIN,
);

if(!$infile){
     print "USAGE: do_saturation_test.pl -i=<INFILE> -bin=<BIN>\n";
     print "<INFILE> input fastq file\n";
     print "<BIN>: number of mapped reads in bin (e.g. 100000)\n";
     exit 1;
}
if(!$INITIAL_BIN){
     print "USAGE: do_saturation_test.pl -i=<INFILE> -bin=<BIN>\n";
     print "<INFILE> input fastq file\n";
     print "<BIN>: number of mapped reads in bin (e.g. 100000)\n";
     exit 1;
}
my $stat = Statistics::Descriptive::Full->new();

#get unique id to customise the random outputs
#TODO this is heavily dependent on the input files and will work only with VDR inputs now
#if($infile =~ /(VDR\_)(\w{2}\d{5})(\_Peconic\d{5})(\_trimmed)(\.\w+)$/){  #eg VDR_GM19213_Peconic20324_trimmed.fastq
#	$input_id = $2;
#}else{
#	print "Input file name not recognised. Modify the regular expression in the script.\n";
#	exit -1;
#}
#TODO REMOVE
$input_id = "vdr_ceu_yri";
#I want all the outputs in a directory whose format will be $PATH/$input_id_saturation_data
my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
if($directory eq "\.\/"){ 
	$datadir = $input_id . '_encode_tests';
}else{
	$datadir =  $directory . "\/" . $input_id . '_encode_tests';
}
system "mkdir $datadir";

#---------
#data - complete
#---------
my $infile_sam              = $datadir . "\/" . $basename . '.sam';
my $infile_bam_predup       = $datadir . "\/" . $basename . '_pd.bam';#this is now kept throughout the whole duration of the script; picard will sample from this
my $infile_bam              = $datadir . "\/" . $basename . '.bam';#this is now kept throughout the whole duration of the script; picard will sample from this
my $infile_bam_dup_metrics  = $datadir . "\/" . $basename . '.bam_metrics';
my $infile_bai              = $datadir . "\/" . $basename . '.bai';
my $infile_bed              = $datadir . "\/" . $basename . '.bed';
my $infile_tagAlign         = $datadir . "\/" . $basename . '.tagAlign';
my $infile_tagAlign_gz      = $datadir . "\/" . $basename . '.tagAlign.gz';
my $infile_tagAlign_pr1     = $datadir . "\/" . $basename . '.pr1.tagAlign.gz';
my $infile_tagAlign_pr2     = $datadir . "\/" . $basename . '.pr2.tagAlign.gz';

my $full_gem_basename       = $datadir . "\/" . $basename  . '_gem';       #gem outputs in these
my $full_macs14_basename    = $datadir . "\/" . $basename  . '_macs14';    #macs14 outputs in these
my $full_macs2_basename     = $datadir . "\/" . $basename  . '_macs2';     #macs2 outputs in these
my $full_macs2_idr_basename = $datadir . "\/" . $basename  . '_macs2idr';  #needed for the macs2 - idr peak call
my $full_macs2_pr1_basename = $datadir . "\/" . $basename  . '_macs2_pr1'; #macs2 outputs here for pseudoreplicate 1
my $full_macs2_pr2_basename = $datadir . "\/" . $basename  . '_macs2_pr2'; #macs2 outputs here for pseudoreplicate 1
my $full_idr_basename       = $datadir . "\/" . $basename  . '_pr1VSpr2';
my $full_idr_plot_basename  = $datadir . "\/" . $basename  . '_IDR_plot';#done only for full, not random

#IDR only, MACS2 only
my $full_macs2_narrowPeak_file        = $full_macs2_idr_basename . '_peaks.narrowPeak';
my $full_macs2_narrowPeak_file_thrs   = $full_macs2_idr_basename . '_peaks-thrs.narrowPeak';
#pseudorep 1 
my $full_macs2_pr1_narrowPeak_file    = $full_macs2_pr1_basename . '_peaks.narrowPeak';
my $full_macs2_pr1_regionPeak_file    = $full_macs2_pr1_basename . '.regionPeak';
my $full_macs2_pr1_regionPeak_gz_file = $full_macs2_pr1_basename . '.regionPeak.gz';
#pseudorep 2
my $full_macs2_pr2_narrowPeak_file    = $full_macs2_pr2_basename . '_peaks.narrowPeak';
my $full_macs2_pr2_regionPeak_file    = $full_macs2_pr2_basename . '.regionPeak';
my $full_macs2_pr2_regionPeak_gz_file = $full_macs2_pr2_basename . '.regionPeak.gz';

my $full_idr_overlap_file             = $full_idr_basename          . '-overlapped-peaks.txt';
my $full_idr_overlap_file_thrs        = $full_idr_basename          . '-overlapped-peaks-thrs.txt';
my $full_idr_encodepeak_file          = $full_idr_overlap_file_thrs . '.npk'; # generated by idrOverlap2npk.sh #input to bedtools multicov
my $full_idr_encodepeak_file_counts   = $full_idr_encodepeak_file   . '.counts'; # needed by do_FRip_IDR
#--------
#data - random
#--------
my $random_bam_file           = $datadir . "\/" . $basename  . '_random.bam'; #picard writes here /macs2 reads
my $random_bai_file           = $datadir . "\/" . $basename  . '_random.bai'; #picard writes here
my $random_bed_file           = $datadir . "\/" . $basename  . '_random.bed'; #bamtobed writes here / genetrack reads
my $random_tagAlign_file      = $datadir . "\/" . $basename  . '_random.tagAlign';
my $random_tagAlign_gz_file   = $datadir . "\/" . $basename  . '_random.tagAlign.gz';
my $random_tagAlign_pr1_file  = $datadir . "\/" . $basename  . '_random.pr1.tagAlign.gz';
my $random_tagAlign_pr2_file  = $datadir . "\/" . $basename  . '_random.pr2.tagAlign.gz';

my $random_gem_basename       = $datadir . "\/" . $basename  . '_random_gem';       #gem outputs in these
my $random_macs14_basename    = $datadir . "\/" . $basename  . '_random_macs14';    #macs14 outputs in these
my $random_macs2_basename     = $datadir . "\/" . $basename  . '_random_macs2';     #macs2 outputs in these
my $random_macs2_idr_basename = $datadir . "\/" . $basename  . '_random_macs2idr';  #needed for the macs2 - idr peak call outputs
my $random_macs2_pr1_basename = $datadir . "\/" . $basename  . '_random_macs2_pr1';  #macs2 outputs here for pseudoreplicate 1
my $random_macs2_pr2_basename = $datadir . "\/" . $basename  . '_random_macs2_pr2';  #macs2 outputs here for pseudoreplicate 1
my $random_idr_basename       = $datadir . "\/" . $basename  . '_random_pr1VSpr2';
my $random_idr_plot_basename  = $datadir . "\/" . $basename  . '_random_IDR_plot';#done only for full, not random

#IDR only, MACS2 only
my $random_macs2_narrowPeak_file        = $random_macs2_idr_basename     . '_peaks.narrowPeak';;
my $random_macs2_narrowPeak_file_thrs   = $random_macs2_idr_basename     . '_peaks-thrs.narrowPeak';
#pseudorep 1
my $random_macs2_pr1_narrowPeak_file    = $random_macs2_pr1_basename . '_peaks.narrowPeak';
my $random_macs2_pr1_regionPeak_file    = $random_macs2_pr1_basename . '.regionPeak';
my $random_macs2_pr1_regionPeak_gz_file = $random_macs2_pr1_basename . '.regionPeak.gz';
#pseudorep 2
my $random_macs2_pr2_narrowPeak_file    = $random_macs2_pr2_basename . '_peaks.narrowPeak';
my $random_macs2_pr2_regionPeak_file    = $random_macs2_pr2_basename . '.regionPeak';
my $random_macs2_pr2_regionPeak_gz_file = $random_macs2_pr2_basename . '.regionPeak.gz';

my $random_idr_overlap_file             = $random_idr_basename          . '-overlapped-peaks.txt';
my $random_idr_overlap_file_thrs        = $random_idr_basename          . '-overlapped-peaks-thrs.txt';
my $random_idr_encodepeak_file          = $random_idr_overlap_file_thrs . '.npk'; # generated by idrOverlap2npk.sh #input to bedtools multicov
my $random_idr_encodepeak_file_counts   = $random_idr_encodepeak_file   . '.counts'; # needed by do_FRip_IDR
my $temp_idr_narrowPeak_file            = $datadir . "\/"               . '_temp.npk';

my $outfile                 = $datadir . "\/" . $basename . '.result';
my $outfile_enrichment_data = $datadir . "\/" . $basename . '_enrichment.result'; # each column in this file is an iteration. Each row is the full signal height data for all new peaks in that iteration. Needed to make boxplots

#unlink files in case something is left from previous run
unlink $infile_sam;
unlink $infile_bam_predup;
unlink $infile_bam;
unlink $infile_bed;
unlink $infile_bai;
unlink $infile_tagAlign_gz;
unlink $infile_tagAlign_pr1;
unlink $infile_tagAlign_pr2;
unlink $random_bam_file;
unlink $random_bed_file;
unlink $random_bai_file;
unlink $random_tagAlign_gz_file;
unlink $random_tagAlign_pr1_file;
unlink $random_tagAlign_pr2_file;

#############################################################
#run pipeline on full dataset to get total number of mapped reads and peaks
#############################################################
my @full_threads; my @full_peaks;	
#------------------
#align (only once)
#------------------
print "-----------------------\n";
print "INFO - Aligning fastq file\n";
print "-----------------------\n";
system "$BOWTIE -q -S --best --strata -m 1 -p 10 --chunkmbs 1024 $INDEX_PATH $infile $infile_sam";
if ( $? == -1 ){
	print "BOWTIE: problem with output: $!\n";
	exit -1;
}
#-------------
#sam to bam
#------------
print "-----------------------\n";
print "INFO - Picard: sam to bam\n";
print "-----------------------\n";
system "java -Xmx4g -Djava.io.tmpdir=$TMP_DIR -jar $PICARD/SortSam.jar SORT_ORDER=coordinate INPUT=$infile_sam OUTPUT=$infile_bam_predup CREATE_INDEX=true VERBOSITY=ERROR VALIDATION_STRINGENCY=LENIENT";
#my $total_mapped_reads = `$SAMTOOLS/samtools view -F 4 $infile_bam_predup | wc -l`; #mapped reads
#chomp $total_mapped_reads;
#-------------
#bam remove duplicates
#------------
print "-----------------------\n";
print "INFO - Picard: bam: remove duplicates\n";
print "-----------------------\n";
system "java -Xmx4g -Djava.io.tmpdir=$TMP_DIR -jar $PICARD/MarkDuplicates.jar INPUT=$infile_bam_predup OUTPUT=$infile_bam METRICS_FILE=$infile_bam_dup_metrics REMOVE_DUPLICATES=true CREATE_INDEX=true VERBOSITY=ERROR VALIDATION_STRINGENCY=LENIENT";
my $total_mapped_reads = `$SAMTOOLS/samtools view -F 4 $infile_bam | wc -l`; #mapped reads
chomp $total_mapped_reads;
#------------------------------
#bam -> tagalign, bam -> bed
#------------------------------
print "-----------------------\n";
print "INFO - bam -> tagalign, bam -> bed\n";
print "-----------------------\n";
my $thread_bam_to_tagalign  = threads->create( \&bam_to_tagalign, $infile_bam);
my $thread_bam_to_bed       = threads->create( \&bam_to_bed, $infile_bam, $infile_bed);
$thread_bam_to_tagalign->join;
$thread_bam_to_bed->join;
if (my $err = $thread_bam_to_tagalign->error()) { warn("Thread bam_to_tagalign() error: $err\n"); } 
if (my $err = $thread_bam_to_bed->error())      { warn("Thread bam_to_bed() error: $err\n");      } 
#-------------------------
#peak calling - unfiltered
#-------------------------
print "-----------------------\n";
print "INFO - Peak Calling - unfiltered\n";
print "-----------------------\n";
my $full_thread_gps       = threads->create( \&run_gem,              $infile_bam,         $full_gem_basename );
my $full_thread_macs14    = threads->create( \&run_macs14,           $infile_tagAlign_gz, $full_macs14_basename );
my $full_thread_macs2     = threads->create( \&run_macs2_get_peaks,  $infile_tagAlign_gz, $full_macs2_basename, $MACS2_OPTS  );
my $full_thread_macs2idr  = threads->create( \&run_macs2_get_peaks,  $infile_tagAlign_gz, $full_macs2_idr_basename, $MACS2_IDR_OPTS  );
#genetrack threads, one per chromosome
for my $i (0 .. $#HUMAN_CHRS) {
	push(@full_threads, threads->create(\&run_genetrack, $HUMAN_CHRS[$i], $infile_bed));
}
foreach my $thread (@full_threads){
	my $peak = $thread->join;
	if (my $err = $thread->error()) { warn("Thread run_genetrack() error: $err\n"); } 
	push(@full_peaks, $peak);	
}
my $total_gem_peaks       = $full_thread_gps->join;
my $total_macs14_peaks    = $full_thread_macs14->join;
my $total_macs2_peaks     = $full_thread_macs2->join;
my $total_macs2idr_peaks  = $full_thread_macs2idr->join; #unused atm
if (my $err = $full_thread_gps->error())    { warn("Thread run_gem() error: $err\n");             }
if (my $err = $full_thread_macs14->error()) { warn("Thread run_macs14() error: $err\n");          } 
if (my $err = $full_thread_macs2->error())  { warn("Thread run_macs2_get_peaks() error: $err\n"); } 
if (my $err = $full_thread_macs2idr->error())  { warn("Thread run_macs2idr_get_peaks() error: $err\n"); } 
my $total_genetrack_peaks = sum(@full_peaks);
#---------------
#FRiP (basic, no IDR)
#---------------
print "-----------------------\n";
print "INFO - get FRiP (GEM, MACS1.4, MACS2)\n";
print "-----------------------\n";
my $total_FRiP_macs14_thread = threads-> create(\&do_FRiP_macs, $total_mapped_reads,  $full_macs14_basename );
my $total_FRiP_macs2_thread  = threads-> create(\&do_FRiP_macs, $total_mapped_reads,  $full_macs2_basename  );
my $total_FRiP_gem_thread    = threads-> create(\&do_FRiP_gem,  $total_mapped_reads,  $full_gem_basename    );
my $total_FRiP_macs14        = $total_FRiP_macs14_thread->join;
my $total_FRiP_macs2         = $total_FRiP_macs2_thread->join;
my $total_FRiP_gem           = $total_FRiP_gem_thread->join;
if (my $err = $total_FRiP_macs14_thread->error())    { warn("Thread frip_macs14() error: $err\n");             }
if (my $err = $total_FRiP_macs2_thread->error()) { warn("Thread frip_macs3() error: $err\n");          } 
if (my $err = $total_FRiP_gem_thread->error())  { warn("Thread frip_gem() error: $err\n"); } 
#------------------------
# IDR: create pseudoreplicates for MACS2 IDR peak calling
#------------------------
print "-----------------------\n";
print "INFO - do_random_split_peakfile.sh: create two random subsets of tagAlign file\n";
print "-----------------------\n";
system "$SCRIPT_PATH/do_random_split_peakfile.sh $infile_tagAlign_gz $datadir";
#-------------------------------
# IDR: call peaks on pseudoreplicates (threaded) with MACS2
# IDR: also call peaks on full file with MACS2 (the peak list will be IDR thresholded)
#-------------------------------
print "-----------------------\n";
print "INFO - IDR - MACS2: Call peaks pseudoreplicates\n";
print "-----------------------\n";
my $full_thread_macs2_pr1  = threads->create( \&run_macs2,  $infile_tagAlign_pr1, $full_macs2_pr1_basename, $MACS2_IDR_OPTS);
my $full_thread_macs2_pr2  = threads->create( \&run_macs2,  $infile_tagAlign_pr2, $full_macs2_pr2_basename, $MACS2_IDR_OPTS);
$full_thread_macs2_pr1->join;
$full_thread_macs2_pr2->join;
if (my $err = $full_thread_macs2_pr1->error()) { warn("Thread run_macs2() error: $err\n"); } 
if (my $err = $full_thread_macs2_pr1->error()) { warn("Thread run_macs2() error: $err\n"); } 
#---------------------------------
#sort MACS2 output narrowPeak file rep1 and rep2 (threaded)
#---------------------------------
print "-----------------------\n";
print "INFO - IDR - sort and gzip MACS2 output narrowPeak file for pseudorep 1 and 2\n";
print "-----------------------\n";
my $full_thread_sort_pr1  = threads->create( \&sort_pseudorep,  $full_macs2_pr1_narrowPeak_file, $full_macs2_pr1_regionPeak_file);
my $full_thread_sort_pr2  = threads->create( \&sort_pseudorep,  $full_macs2_pr2_narrowPeak_file, $full_macs2_pr2_regionPeak_file);
$full_thread_sort_pr1->join;
$full_thread_sort_pr2->join;
if (my $err = $full_thread_sort_pr1->error()) { warn("Thread sort_pseudorep() error: $err\n"); } 
if (my $err = $full_thread_sort_pr2->error()) { warn("Thread sort_pseudorep() error: $err\n"); } 
#------------
#IDR Analysis
#------------
print "-----------------------\n";
print "INFO - IDR run on pseudoreplicates\n";
print "-----------------------\n";
system "$R_CODE/Rscript $IDR/batch-consistency-analysis.r $full_macs2_pr1_regionPeak_file $full_macs2_pr2_regionPeak_file -1 $full_idr_basename 0 F p.value";
#consistency plot (done only for the full file, not the random subsets later)
system "$R_CODE/Rscript $IDR/batch-consistency-plot.r 1 $full_idr_plot_basename  $full_idr_basename";
my $total_IDR_peaks = `$SCRIPT_PATH/idr_get_peaks.sh  $full_idr_overlap_file $IDR_THRESHOLD | wc -l`; #get the number of peak which pass the threshold
chomp $total_IDR_peaks;
#---
print "Genetrack: number of total peaks in initial file: "           . $total_genetrack_peaks  . "\n";
print "GPS: number of total peaks in initial file: "                 . $total_gem_peaks        . "\n";
print "MACS14: number of total peaks in initial file: "              . $total_macs14_peaks     . "\n";
print "MACS2: number of total peaks in initial file: "               . $total_macs2_peaks      . "\n";
print "MACS2-IDR: number of total filtered peaks in initial file: "  . $total_IDR_peaks        . "\n";
#----------------
#IDR - get candidate peak list
#----------------
system "cat $full_macs2_narrowPeak_file | sort -k8nr,8nr | head -n $total_IDR_peaks > $full_macs2_narrowPeak_file_thrs";


#get subset of -overlapped-peaks.txt that pass $IDR_THRESHOLD ($full_idr_overlap_file_thrs)
#convert $full_idr_overlap_file_thrs output to .narrowPeak format
#use slightly modified version of anshul's script
#USAGE: idrOverlap2npk.sh [idrOverlapFile] [oDir]
#odir optional
#----------------------
#system "$SCRIPT_PATH/idr_get_peaks.sh $full_idr_overlap_file $IDR_THRESHOLD > $full_idr_overlap_file_thrs"; #get the actual peaks which pass the threshold
#system "$SCRIPT_PATH/idrOverlap2npk.sh $full_idr_overlap_file_thrs";
#output file is $full_idr_encodepeak_file

#in order to obtain the FRiP from the IDR peaks I need to count the reads in peak first - can be done with bedtools
#bedtools multicov -bams aln1.bam -bed ivls-of-interest.bed
#the narrowpeak file converted from idr output via idrOverlap2npk SHOULD work (otherwise you need to get a bed directly)
#system "$BEDTOOLS/multiBamCov -bams $infile_bam -bed $full_idr_encodepeak_file  > $full_idr_encodepeak_file_counts"; this works only if you rename the .bai file to .bam.bai
#system "$BEDTOOLS/bedtools coverage -abam $infile_bam -b $full_idr_encodepeak_file -counts > $full_idr_encodepeak_file_counts";
#extending the intervals left and right by SLOP bp (to get more reads in the peaks)
system "$BEDTOOLS/bedtools slop -g $CHROM_SIZE_PATH -b $SLOP -i $full_macs2_narrowPeak_file_thrs | sort -k1,1V -k2,2g | $BEDTOOLS/bedtools merge -i stdin | $BEDTOOLS/bedtools coverage -abam $infile_bam -b stdin -counts > $full_idr_encodepeak_file_counts";
#output is now in format
#chr1    564445  564721  138
#chr1    569806  570063  167
#chr1    2068622 2068748 2
#chr1    3713073 3713221 4
#...
my $total_FRiP_IDR    = do_FRiP_IDR ($total_mapped_reads, $full_idr_encodepeak_file_counts);

#remove all temp files apart from the bams (full/ noduplicates)
unlink $infile_sam; 
unlink $infile_bai;
unlink $infile_bed;
unlink $infile_tagAlign_gz;
unlink $infile_tagAlign_pr1;
unlink $infile_tagAlign_pr2;
my $gt_outputs_gff = $datadir . "\/" . "$input_id\*chr\*.gff";
my $gt_outputs_bed = $datadir . "\/" . "$input_id\*chr\*.bed";
system "rm $gt_outputs_gff"; #rm genetrack outputs
system "rm $gt_outputs_bed"; #rm genetrack outputs
system "rm -rf $full_gem_basename*"; #remove gem output
system "rm -rf $full_macs14_basename*"; #remove macs14 output
system "rm -rf $full_macs2_basename*"; #remove macs2 output
system "rm -rf $full_idr_basename*"; #remove IDR output (but not the plot)

#########################
#randomisation main loop
#########################
#outstream header
open (my $outstream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
open (my $outstream_enrichment,  q{>}, $outfile_enrichment_data) or die("Unable to open $outfile_enrichment_data : $!");
print $outstream "S\tM_R\tP_GENETRACK\tPR_GENETRACK\tP_GEM\tPR_GEM\tP_M14\tPR_14\tP_M2\tPR_M2\tP_IDR\tPR_IDR\tFRIP_M14\tFRIP_M2\tFRIP_GEM\tFRIP_IDR\tM_E\n";
#print $outstream "SAMPLE\tMAPPED_READS\tPEAKS_GENETRACK\tPEAK_RATIO_GENETRACK\tPEAKS_GPS\tPEAK_RATIO_GPS\tPEAKS_MACS14\tPEAK_RATIO_MACS14\tPEAKS_MACS2\tPEAK_RATIO_MACS2\tPEAKS_IDR_MACS2\tPEAK_RATIO_IDR_MACS2\tFRIP_GPS\tFRIP_IDR\tMEDIAN_ENRICHMENT\n";
system "touch $temp_idr_narrowPeak_file";#initally empty, do control first time
my $counter = 1;
my $BIN = $INITIAL_BIN;
while ($BIN <= $total_mapped_reads){
	my $wc_out; my @wc_out;
	my @gt_threads; my @gt_peaks;
	#----------------
	#random subset selection from bam file (Picard)
	#----------------
	my $p_rounded = sprintf("%.3f", ( $BIN / $total_mapped_reads ) );
	print "Creating random .bam file, (fraction:  ( $BIN / $total_mapped_reads ) )..\n";
    system "java -Xmx4g -Djava.io.tmpdir=$TMP_DIR -jar $PICARD/DownsampleSam.jar INPUT=$infile_bam OUTPUT=$random_bam_file RANDOM_SEED=null VERBOSITY=ERROR PROBABILITY=$p_rounded";
    my $number_of_mapped_reads = `$SAMTOOLS/samtools view -F 4 $random_bam_file | wc -l`; #uniquely mapped reads
	chomp $number_of_mapped_reads;
	#------------------------------
	#bam -> tagalign and bam -> bed (threaded)
	#------------------------------
	print "-----------------------\n";
	print "INFO - bam -> tagalign, bam -> bed\n";
	print "-----------------------\n";
	my $random_thread_bam_to_tagalign  = threads->create( \&bam_to_tagalign, $random_bam_file);
	my $random_thread_bam_to_bed       = threads->create( \&bam_to_bed,      $random_bam_file, $random_bed_file);
	$random_thread_bam_to_tagalign->join;
	$random_thread_bam_to_bed->join;
	if (my $err = $random_thread_bam_to_tagalign->error()) { warn("Thread bam_to_tagalign() error: $err\n"); } 
	if (my $err = $random_thread_bam_to_bed->error()) { warn("Thread bam_to_bed() error: $err\n"); } 
	#-------------------------
	#peak calling - unfiltered
	#-------------------------
	print "-----------------------\n";
	print "INFO - Peak Calling - unfiltered\n";
	print "-----------------------\n";
	my $tgps       = threads->create( \&run_gem,              $random_bam_file, $random_gem_basename);
	my $tmacs14    = threads->create( \&run_macs14,           $random_tagAlign_gz_file, $random_macs14_basename );
	my $tmacs2     = threads->create( \&run_macs2_get_peaks,  $random_tagAlign_gz_file, $random_macs2_basename,     $MACS2_OPTS  );
	my $tmacs2idr  = threads->create( \&run_macs2_get_peaks,  $random_tagAlign_gz_file, $random_macs2_idr_basename, $MACS2_IDR_OPTS  );
	#genetrack threads, one per chromosome
	for my $i (0 .. $#HUMAN_CHRS) {
		push(@gt_threads, threads->create(\&run_genetrack, $HUMAN_CHRS[$i], $random_bed_file));
	}
	#push (@gt_peaks, $_->join) foreach @gt_threads;
	foreach my $thread (@gt_threads){
    	my $peak = $thread->join;
		if (my $err = $thread->error()) { warn("Thread run_genetrack() error: $err\n"); } 
		push(@gt_peaks, $peak);	
	}
	my $macs14_peaks = $tmacs14->join;
	my $macs2_peaks  = $tmacs2->join;
	my $macs2idr_peaks = $tmacs2idr->join; #this output is unused currently (you use other macs2 output)
	my $gem_peaks    = $tgps->join;
	if (my $err = $tmacs14->error()) { warn("Thread run_macs14() error: $err\n"); }
	if (my $err = $tmacs2->error()) { warn("Thread run_macs2_get_peaks() error: $err\n"); }
	if (my $err = $tmacs2idr->error()) { warn("Thread run_macs2_get_peaks() error: $err\n"); }
	if (my $err = $tgps->error()) { warn("Thread run_gem() error: $err\n"); }
	my $genetrack_peaks = sum(@gt_peaks);
	#------------------------
	# create pseudoreplicates
	#------------------------
	print "-----------------------\n";
	print "INFO - do_random_split_peakfile.sh: create two random subsets of random tagAlign sample\n";
	print "-----------------------\n";
	system "$SCRIPT_PATH/do_random_split_peakfile.sh $random_tagAlign_gz_file $datadir";
	#-------------------------------
	# call peaks on pseudoreplicates (threaded)
	#-------------------------------
	print "-----------------------\n";
	print "INFO - MACS2: Call peaks pseudoreplicates\n";
	print "-----------------------\n";
	my $random_thread_macs2_pr1  = threads->create( \&run_macs2,  $random_tagAlign_pr1_file, $random_macs2_pr1_basename, $MACS2_IDR_OPTS);
	my $random_thread_macs2_pr2  = threads->create( \&run_macs2,  $random_tagAlign_pr2_file, $random_macs2_pr2_basename, $MACS2_IDR_OPTS);
	$random_thread_macs2_pr1->join;
	$random_thread_macs2_pr2->join;
	if (my $err = $random_thread_macs2_pr1->error()) { warn("Thread run_macs2() error: $err\n"); }
	if (my $err = $random_thread_macs2_pr2->error()) { warn("Thread run_macs2() error: $err\n"); }
	#---------------------------------
	#sort MACS2 output narrowPeak file rep1 and rep2 (threaded)
	#---------------------------------
	print "-----------------------\n";
	print "INFO - sort and gzip MACS2 output regionPeak file for pseudorep 1 and 2\n";
	print "-----------------------\n";
	my $random_thread_sort_pr1  = threads->create( \&sort_pseudorep,  $random_macs2_pr1_narrowPeak_file, $random_macs2_pr1_regionPeak_file);
	my $random_thread_sort_pr2  = threads->create( \&sort_pseudorep,  $random_macs2_pr2_narrowPeak_file, $random_macs2_pr2_regionPeak_file);
	$random_thread_sort_pr1->join;
	$random_thread_sort_pr2->join;
	if (my $err = $random_thread_sort_pr1->error()) { warn("Thread sort_pseudorep() error: $err\n"); }
	if (my $err = $random_thread_sort_pr2->error()) { warn("Thread sort_pseudorep() error: $err\n"); }
	#-------------
	#FRiP - basic
	#-------------	
	print "-----------------------\n";
	print "INFO - FRiP (Macs 1.4, Macs 2, GPS)\n";
	print "-----------------------\n";
	my $random_FRiP_macs14 = do_FRiP_macs($number_of_mapped_reads, $random_macs14_basename);
	my $random_FRiP_macs2  = do_FRiP_macs($number_of_mapped_reads, $random_macs2_basename);
	my $random_FRiP_gem    = do_FRiP_gem ($number_of_mapped_reads, $random_gem_basename);
	#------------
	#IDR Analysis
	#------------
	print "-----------------------\n";
	print "INFO - IDR run on pseudoreplicates\n";
	print "-----------------------\n";
	system "$R_CODE/Rscript $IDR/batch-consistency-analysis.r $random_macs2_pr1_regionPeak_file $random_macs2_pr2_regionPeak_file -1 $random_idr_basename 0 F p.value";
	my $IDR_peaks = `$SCRIPT_PATH/idr_get_peaks.sh  $random_idr_overlap_file $IDR_THRESHOLD | wc -l`;
	chomp $IDR_peaks;
	#------------------------
	#IDR - get candidate peaks
	#------------------------
	system "cat $random_macs2_narrowPeak_file | sort -k8nr,8nr | head -n $IDR_peaks > $random_macs2_narrowPeak_file_thrs";
	
	#system "$SCRIPT_PATH/idr_get_peaks.sh  $random_idr_overlap_file $IDR_THRESHOLD > $random_idr_overlap_file_thrs";
	#system "$SCRIPT_PATH/idrOverlap2npk.sh $random_idr_overlap_file_thrs"; #output is in $random_idr_encodepeak_file
	#----------------------
	#IDR - median fold enrichment (using only IDR-thresholded output peaks) 
	#(since there is no control this is the median value in field 7 of the .narrowPeak file)
	#----------------------
	#my $median_newpeak_signal = get_median_peaksignal_for_newpeaks($random_idr_encodepeak_file, $temp_idr_narrowPeak_file);
	my $median_newpeak_signal = get_median_peaksignal_for_newpeaks($random_macs2_narrowPeak_file_thrs, $temp_idr_narrowPeak_file);
	system "rm $temp_idr_narrowPeak_file"; 
	#system "cp $random_idr_encodepeak_file $temp_idr_narrowPeak_file"; #copy encode peak in temp file for next iteration
	system "cp $random_macs2_narrowPeak_file_thrs $temp_idr_narrowPeak_file"; #copy encode peak in temp file for next iteration
	print "MACS2: number of unfiltered peaks in sample "       . $counter . ": "  . $macs2_peaks           . "\n";
	print "IDR: number of total IDR filtered peaks in sample " . $counter . ": "  . $IDR_peaks             . "\n";
	print "Median New peak signal "                            . $counter . ": "  . $median_newpeak_signal . "\n";
	#----------------------
	#IDR - FRiP
	#----------------------	
	#in order to obtain the FRiP from the IDR peaks I need to count the reads in peak first
	#the narrowpeak file converted from idr output via idrOverlap2npk SHOULD work (otherwise you need to get a bed directly)
 	#system "$BEDTOOLS/multiBamCov -bams $random_bam_file -bed $random_idr_encodepeak_file  > $random_idr_encodepeak_file_counts"; #this works only when the bam index has the name in samtools format, but I use picard :(
 	#TODO: BENEFICIAL? extend intervals by X bp left and right
 	#bedtools slop [OPTIONS] -i <BED/GFF/VCF> -g <GENOME> [-b or (-l and -r)]
 	#system "$BEDTOOLS/bedtools coverage -abam $random_bam_file -b $random_idr_encodepeak_file -counts > $random_idr_encodepeak_file_counts";
	system "$BEDTOOLS/bedtools slop -g $CHROM_SIZE_PATH -b $SLOP -i $random_macs2_narrowPeak_file_thrs | sort -k1,1V -k2,2g | $BEDTOOLS/bedtools merge -i stdin | $BEDTOOLS/bedtools coverage -abam $random_bam_file -b stdin -counts > $random_idr_encodepeak_file_counts";
	#TODO check the output should be the same file plus a row with the numbers. Sum those numbers
	my $random_FRiP_IDR    = do_FRiP_IDR ($number_of_mapped_reads, $random_idr_encodepeak_file_counts);
	
	print $outstream  $counter .
	                  "\t" . $number_of_mapped_reads .
	                  "\t" . $genetrack_peaks  . "\t" . ($genetrack_peaks/$total_genetrack_peaks) .  
	                  "\t" . $gem_peaks        . "\t" . ($gem_peaks/$total_gem_peaks) .
	                  "\t" . $macs14_peaks     . "\t" . ($macs14_peaks/$total_macs14_peaks) .
					  "\t" . $macs2_peaks      . "\t" . ($macs2_peaks/$total_macs2_peaks) .
					  "\t" . $IDR_peaks        . "\t" . ($IDR_peaks/$total_IDR_peaks) .
					  "\t" . $random_FRiP_macs14    .
					  "\t" . $random_FRiP_macs2     .
					  "\t" . $random_FRiP_gem       .
					  "\t" . $random_FRiP_IDR       .
					  "\t" . $median_newpeak_signal .
					  "\n" ;

	$BIN += $INITIAL_BIN;
	$counter += 1;
	unlink $random_bam_file;
	unlink $random_bai_file;
	unlink $random_bed_file;
	unlink $random_tagAlign_gz_file;
	unlink $random_tagAlign_pr1_file;
	unlink $random_tagAlign_pr2_file;

	system "rm $gt_outputs_gff"; #rm genetrack outputs
	system "rm $gt_outputs_bed"; #rm genetrack outputs
    system "rm $random_macs14_basename*.*"; #rm macs14 output
	system "rm $random_macs2_basename*.*"; #rm macs2 output
	system "rm -rf $random_gem_basename*";#rm gem output
	system "rm -rf $random_idr_basename*"; #remove IDR output
}
print $outstream  $counter 
                  . "\t" . $total_mapped_reads
                  . "\t" . $total_genetrack_peaks  . "\t" . ($total_genetrack_peaks/$total_genetrack_peaks)  
                  . "\t" . $total_gem_peaks        . "\t" . ($total_gem_peaks/$total_gem_peaks) 
                  . "\t" . $total_macs14_peaks     . "\t" . ($total_macs14_peaks/$total_macs14_peaks)
				  . "\t" . $total_macs2_peaks      . "\t" . ($total_macs2_peaks/$total_macs2_peaks) 
				  . "\t" . $total_IDR_peaks        . "\t" . ($total_IDR_peaks/$total_IDR_peaks)
		          . "\t" . $total_FRiP_macs14
				  . "\t" . $total_FRiP_macs2
				  . "\t" . $total_FRiP_gem
				  . "\t" . $total_FRiP_IDR
				  . "\t" . 'NA'
				  . "\n" ;
unlink $infile_bam;
unlink $infile_bam_predup;
system "rm $temp_idr_narrowPeak_file";
close $outstream;
close $outstream_enrichment;
print "\nFINISHED\n";

##################
# SUBROUTINES
##################
#-----------------------------------------------------------------------------
#bedtools subtract -a genes.bed -b introns.bed
#the -A option will completely remove a feature from A if it has even 1bp of overlap with a feature in B
#warning - maybe it doesn't work with encodepeaks files
#chr11	0	134946492	GM00000_saturation_test/VDR_GM00000_Peconic00000_trimmed_macs2_peak_2	76	.	1.99957	7.63337	-1.00000	134946479
sub get_median_peaksignal_for_newpeaks{
	my  ( $narrowPeak_file, $previous_narrowPeak_file ) = @_;
	my $working_peak_file = $datadir . "\/" .  '_working.npk'; ;
	#system "rm $working_peak_file";
	my @novel_peak_signals;

	if (-s $previous_narrowPeak_file) { #if the temp file is not empty
		system "$BEDTOOLS/subtractBed -A -a $narrowPeak_file -b $previous_narrowPeak_file  > $working_peak_file";
	}else{
		system "cp $narrowPeak_file $working_peak_file";#it means this is the first iteration	
	}
	open (my $peak_stream,       q{<}, $working_peak_file) or die ("Unable to open $working_peak_file: $!");
	while(<$peak_stream>){
        chomp($_);
		my $peak_signal = (split /\t/)[6];
		#print $peak_signal, "\n";
		push(@novel_peak_signals, $peak_signal);
	}
	close $peak_stream;
	system "rm $working_peak_file";

	#create row in output file of total new peak signals
	print $outstream_enrichment join("\t", @novel_peak_signals), "\n";
	
	#get only median for main output file
	$stat->add_data(@novel_peak_signals);
	return $stat->median();
}
	
#-----------------------------------------------------------------------------
sub bam_to_tagalign{
	my ($bam_input) = @_;
	system "$SCRIPT_PATH/do_bam_to_tagalign.sh $bam_input $datadir";
	if ( $? == -1 ){
    	print "bam_to_tagalign(): problem with output: $!\n";
    	exit -1;
	}
	return;
}
#-----------------------------------------------------------------------------
sub bam_to_bed{
	my ($bam_input, $bed_output) = @_;
	system "$BEDTOOLS/bamToBed -i $bam_input > $bed_output";
	if ( $? == -1 ){
    	print "bam_to_bed(): problem with output: $!\n";
    	exit -1;
	}
	return;
}
#-----------------------------------------------------------------------------
sub sort_pseudorep{
	my ( $narrowPeak_file, $regionPeak_file ) = @_;
	system "sort -k 8nr,8nr $narrowPeak_file | head -n 100000 > $regionPeak_file";
	return;
}
#-----------------------------------------------------------------------------
#-------------
# PEAK CALLERS
#-------------
#genetrack output is as follows
#chr1    genetrack       .       568373  568383  6.24515535516   +       .       readcount=6;ID=P568378;stddev=2.35702260396;height=6.24515535516
#chr1    genetrack       .       568477  568487  6.0292452985    +       .       readcount=5;ID=P568482;stddev=0.0;height=6.0292452985
sub run_genetrack{
	my ( $chr, $input_file ) = @_; 
	my $temp_file        = $datadir . "\/" . $input_id .  '_'  . $chr . '_' . 'raw.gff';
	my $gff_file         = $datadir . "\/" . $input_id .  '_'  . $chr .       '.gff';
	my $sorted_gff_file  = $datadir . "\/" . $input_id .  '_'  . $chr . '_' . 'sorted.gff';
	my $sorted_bed_file  = $datadir . "\/" . $input_id .  '_'  . $chr .       '.bed'; # this is the one where you need to count

	#python ${PCODE}/genetrack.py -v -s 5 -e 10 -F 6 -c chr1 ${FILE} #fine grain
	#system "python $GENETRACK -s 10 -e 40 -c $chr $input_file > $temp_file"; #course grain
	#system "source activate";
	system "python $GENETRACK -s 5 -e 10 -c $chr $input_file > $temp_file";
    if ( $? == -1 ){
    	print "Genetrack: problem with output: $!\n";
        exit -1;
    }
	#1 - remove all "singletons" (entries where sd = 0.0)
	open (my $temp_data,  q{<}, $temp_file) or die ("Unable to open $temp_file: $!");
	open (my $out_data,  q{>}, $gff_file) or die ("Unable to open file $gff_file: $!");
	while(<$temp_data>){
		print $out_data $_ unless($_ =~ /stddev=0.0;/);
	}
	close $temp_data;
	close $out_data;
	#2 - merge + and - peaks
	#sort by chr and by start coordinate first (required by bedtools)
	system "sort -k1,1 -k4n,4 $gff_file > $sorted_gff_file";
	if ( $? == -1 ){
		print "sort: problem with output: $!\n";
		exit -1;
	}
	#YOU NEED TO OPTIMISE THAT -d!!
	system "$BEDTOOLS/mergeBed -n -d 35 -i $sorted_gff_file > $sorted_bed_file";
	if ( $? == -1 ){
		print "mergeBed: problem with output: $!\n";
		exit -1;
	}
	my $wc_out = `wc -l $sorted_bed_file`;
    my @peaks = split(/ /, $wc_out);
	#$peaks[0] -= 1 if ($peaks[0] > 0); #genetrack produces 1 line of header
	return $peaks[0];
}
#-----------------------------------------------------------------------------
sub run_gem{
	my ( $input_file, $gem_basename ) = @_;
	system "java -Xmx10G -jar $GEM/gem.jar --g $CHROM_SIZE_PATH --d $GEM/Read_Distribution_ChIP-exo.txt --s $MAPPABILITY --expt $input_file --f SAM  --out $gem_basename $GEM_OPTS > /dev/null";
	#system "java -Xmx10G -jar $GEM/gem.jar --g $CHROM_SIZE_PATH --d $GEM/Read_Distribution_ChIP-exo.txt --s $MAPPABILITY --expt $input_file --f SAM --top 0 --out $gem_basename";
	if ( $? == -1 ){
                print "GEM: problem with output: $!\n";
                exit -1;
    }
    #now capture number of rows in bedfile and fill an output file with them
    #TODO change this if you use GPS instead of gem
    my $gem_peakfile = $gem_basename . "_GPS_events.txt";
    #my $gem_peakfile = $gem_basename . "_GEM_events.txt";
    my $wc_out = `wc -l $gem_peakfile`;
    my @gem_peaks = split(/ /, $wc_out);
	return $gem_peaks[0];
}
#-----------------------------------------------------------------------------
sub run_macs14{
	my ( $input_file, $macs14_basename ) = @_;
	#/macs14 -t ${FILE} -n ${PDATA}/${ID} -s 40 --nomodel --shiftsize=13 --keep-dup=all -p 1e-8 -g 2540757438
	#system "macs14 -t  $input_file -n $macs14_basename -g $MAPPABILITY -s 40 --keep-dup=all --nomodel --shiftsize 13 --verbose 1 --nolambda --llocal 0";
	#system "source activate";
	system "macs14 -t  $input_file -f BED -n $macs14_basename -g $MAPPABILITY $MACS14_OPTS";	
	if ( $? == -1 ){
		print "MACS14: problem with output: $!\n";
                exit -1;
        }
	#now capture number of rows in bedfile and fill an outputfile with the
	my $macs14_bedfile = $macs14_basename . "_peaks.bed";
    my $wc_out = `wc -l $macs14_bedfile`;
    my @peaks = split(/ /, $wc_out);
	return $peaks[0];	
}
#-----------------------------------------------------------------------------
sub run_macs2_get_peaks{
	my ( $input_file, $macs2_basename, $opts ) = @_;
	
	run_macs2($input_file, $macs2_basename, $opts);
	#now capture number of rows in narrowPeak bed file
	my $macs2_bedfile = $macs2_basename . "_peaks.narrowPeak";
    my $wc_out = `wc -l $macs2_bedfile`;
    my @peaks = split(/ /, $wc_out);
	return $peaks[0];
}
#-----------------------------------------------------------------------------
sub run_macs2{
	my ( $input_file, $macs2_basename, $opts ) = @_;
	#${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/${ID} -g 2540757438 --keep-dup=all -q 0.001 --nomodel --extsize 26 --call-summits
	#system "macs2 callpeak -t $input_file -n $macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 --verbose 1 --nolambda --llocal 0";
	#use the following from IDR:
	#macs2 callpeak -t $infile_tagAlign_gz -f BED -n $full_macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 -p 1e-3 --to-large
	#system "macs2 callpeak -t $input_file -f BED -n $macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 -p 1e-3 --verbose 1 --nolambda --llocal 0";
	#system "macs2 callpeak -t $input_file -f BED -n $macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 -p 1e-3 --verbose 1";
	#system "source activate";
	system "macs2 callpeak -t $input_file -f BED -n $macs2_basename -g $MAPPABILITY $opts";
	if ( $? == -1 ){
		print "MACS2: problem with output: $!\n";
        exit -1;
    }
    return;	
}
#-----------------------------------------------------------------------------
##################
#FRiP
##################
#-----------------------------------------------------------------------------
sub do_FRiP_gem{
	my ( $total_reads, $basename ) = @_;
	#open even file
	#TODO change this if you use GPS rather than GEM
	#my $peak_data = $basename . "_GPS_events.txt";
	my $peak_data = $basename . "_GEM_events.txt";
	my @pileup;	my $total_reads_in_peaks;
	
	open (my $peak_stream,  q{<}, $peak_data) or die ("Unable to open $peak_data: $!");
	while(<$peak_stream>){
		chomp;
		next if($_ =~ /^Position/); #header	
			
		if($_ =~ /^\w+\:\d+/){  #1:565745 #riga dati
			if($_ =~ /^(\w+\:\d+)\s+(\d+\.\d+)\s+(NaN)/){ #1:565745          298.4 NaN
				my $read_count = $2;
				push(@pileup, $read_count);
			}else{
				print "it shouldn't be here, cannot use regex\n";
				exit -1;
			}
		}else{
			print ("do_FRiP_gem(): error, unrecognised line in .txt: $_\n");
			exit -1;
		}
	}
	close $peak_stream;
	$total_reads_in_peaks = sum(@pileup);
	return ($total_reads_in_peaks/$total_reads);
}
#-----------------------------------------------------------------------------
#this gets an .xls file produced my macs, and the total number or mapped reads
#for eack peak in the xls, it gets the "pileup" column, sums for all peaks
#it returns sum(pileup)/(total number of aligned reads)
sub do_FRiP_macs{
	my ( $total_reads, $basename ) = @_;
	#open excel file
	my $peak_data = $basename . '_peaks.xls'; 
	my @pileup;	my $total_reads_in_peaks;
	
	open (my $peak_stream,  q{<}, $peak_data) or die ("Unable to open $peak_data: $!");
	while(<$peak_stream>){
		chomp;
		next if($_ =~ /^\#/); #macs comment
		next if($_ eq ''); #blank line
		next if($_ =~ /chr\t/); #header	
		
		if($_ =~ /^chr\w+/){
			my $read_count = (split /\t/)[5];
			push(@pileup, $read_count);
		}else{
			print ("do_FRiP(): error, unrecognised line in .xls: $_\n");
			exit -1;
		}
	}
	close $peak_stream;
	$total_reads_in_peaks = sum(@pileup);
	#print 'MACS: Total reads in peaks: ', $total_reads_in_peaks, "\n";
	return ($total_reads_in_peaks/$total_reads);
}
#-----------------------------------------------------------------------------
#peak file format (npk + bedtools count)
#chr1    0       565729  1       2311.39898      .       2083.05609      4.426070        4.426070        1900
sub do_FRiP_IDR{
	my ( $total_reads, $peak_data ) = @_;
	my @pileup;
	my $total_reads_in_peaks;
	
	open (my $peak_stream,  q{<}, $peak_data) or die ("Unable to open $peak_data: $!");
	while(<$peak_stream>){
		chomp;
		my $read_count = (split /\t/)[3];
	#	print $read_count, "\n";
		push(@pileup, $read_count);
	}
	close $peak_stream;
	$total_reads_in_peaks = sum(@pileup);
	return ($total_reads_in_peaks/$total_reads);
}
