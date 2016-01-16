#!/usr/bin/perl
use strict;
use warnings;
use threads;
use Getopt::Long;
use File::Basename;
use Statistics::Descriptive;
use List::Util qw(sum);

my $stat = Statistics::Descriptive::Full->new();

#13/3/2013
#This is obtained from saturation_test_mt_4.pl and marginal_fe_test.pl
#1 does both saturation_analysis and median enrichment plots
#2 for the saturation analysis, it does genetrack(+postproc), GPS, MACS 1.4, MACS2 (+IDR)
#does the full postprocessing IDR pipeline (standalone in the _idr_pseudorep.pl file) after each MACS2 peak calling

#IDR is very VERY slow and the other threads will all be waiting for the idr ones
#also Anshul suggests not to use on MACS14
#should probably also try it with GEM
#should probably parallelise the random subsets, for now they happen in sequence and all output for one iteration is deleted after output peak number is collected

#INPUTS
#1 fastq file
#2 step or frequency (0.01? 0.1)

#The step approach is wrong, because I will get less samples the larger the file.
#If I remove GM10846, GM12752 and GM11919 from the computation (all less than 1.000.000 uniquely mapped reads) I have read numbers ranging from  1.000.000 to 15.000.000
#I could try one iteration every 200.000 reads

#OUTPUTS
#TSV FILE  <samplenumber> <uniquely mapped reads> [<number of peaks> <number of peaks/total peaks for full file>] <marginal_fe>

my $infile;my $INITIAL_BIN;
my $input_id; my $datadir; # to customise randomised outputs
my $SCRIPT_PATH = '/net/isi-scratch/giuseppe/scripts';
my $TOOL_PATH = '/net/isi-scratch/giuseppe/tools';
my $INDEX_PATH = '/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19_bowtie';
my $CHROM_SIZE_PATH = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes';
my $GENOME_PATH = '/net/isi-scratch/giuseppe/indexes/Hsap/hg19_genome';
my $R_CODE = '/net/isi-cgat/ifs/apps/apps/R-2.14.1/bin';
my $TMP_DIR = '/tmp';

#program paths, change as needed
my $BOWTIE    = $TOOL_PATH . '/bowtie-0.12.9/bowtie';
my $PICARD    = $TOOL_PATH . '/picard-tools-1.84/picard-tools-1.84';
my $SAMTOOLS  = $TOOL_PATH . '/samtools';
my $BEDTOOLS  = $TOOL_PATH . '/bedtools-2.17.0/bin';
my $IDR       = $TOOL_PATH . '/idrCode'; #now using idr in giuseppe/scripts
my $GEM       = $TOOL_PATH . '/gem';
my $GENETRACK = $TOOL_PATH . '/chipexo-master/genetrack/genetrack.py';
#my $SPP = $TOOL_PATH . '/phantompeakqualtools';

#human- and dataset-specific
my $MAPPABILITY = 2540757438; #for Hsap (40bp) used by Macs2 and GEM
my @HUMAN_CHRS = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'); 

GetOptions(
	'input=s'  =>\$infile,
	'bin=f' => \$INITIAL_BIN,
);

if(!$infile){
     print "USAGE: do_saturation_test.pl -input=<INFILE> -bin=<BIN>\n";
     print "<INFILE> input fastq file\n";
     print "<BIN>: number of mapped reads in bin (e.g. 100000)\n";
     exit 1;
}
if(!$INITIAL_BIN){
     print "USAGE: do_saturation_test.pl -input=<INFILE> -bin=<BIN>\n";
     print "<INFILE> input fastq file\n";
     print "<BIN>: number of mapped reads in bin (e.g. 100000)\n";
     exit 1;
}
#get unique id to customise the random outputs
#TODO this is heavily dependent on the input files and will work only with VDR inputs now
if($infile =~ /(VDR\_)(\w{2}\d{5})(\_Peconic\d{5})(\_trimmed)(\.\w+)$/){  #eg VDR_GM19213_Peconic20324_trimmed.fastq
	$input_id = $2;
}else{
	print "Input file name not recognised. Modify the regular expression in the script.\n";
	exit -1;
}
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
my $infile_bam              = $datadir . "\/" . $basename . '.bam';#this is now kept throughout the whole duration of the script; picard will sample from this
my $infile_bai              = $datadir . "\/" . $basename . '.bai';
my $infile_bed              = $datadir . "\/" . $basename . '.bed';
my $infile_tagAlign         = $datadir . "\/" . $basename . '.tagAlign';
my $infile_tagAlign_gz      = $datadir . "\/" . $basename . '.tagAlign.gz';
my $infile_tagAlign_pr1     = $datadir . "\/" . $basename . '.pr1.tagAlign.gz';
my $infile_tagAlign_pr2     = $datadir . "\/" . $basename . '.pr2.tagAlign.gz';

my $full_gem_basename       = $datadir . "\/" . $basename  . '_gem';       #gem outputs in these
my $full_macs14_basename    = $datadir . "\/" . $basename  . '_macs14';    #macs14 outputs in these
my $full_macs2_basename     = $datadir . "\/" . $basename  . '_macs2';      #macs2 outputs in these
my $full_macs2_pr1_basename = $datadir . "\/" . $basename  . '_macs2_pr1';  #macs2 outputs here for pseudoreplicate 1
my $full_macs2_pr2_basename = $datadir . "\/" . $basename  . '_macs2_pr2';  #macs2 outputs here for pseudoreplicate 1
my $full_idr_basename       = $datadir . "\/" . $basename  . '_pr1VSpr2';
my $full_idr_plot_basename  = $datadir . "\/" . $basename  . '_IDR_plot';#done only for full, not random

#IDR only, MACS2 only
my $full_macs2_encodePeak_file        = $full_macs2_basename . '_peaks.encodePeak';
my $full_macs2_regionPeak_file        = $full_macs2_basename . '.regionPeak';
my $full_macs2_regionPeak_gz_file     = $full_macs2_basename . '.regionPeak.gz';
#pseudorep 1 
my $full_macs2_pr1_encodePeak_file    = $full_macs2_pr1_basename . '_peaks.encodePeak';
my $full_macs2_pr1_regionPeak_file    = $full_macs2_pr1_basename . '.regionPeak';
my $full_macs2_pr1_regionPeak_gz_file = $full_macs2_pr1_basename . '.regionPeak.gz';
#pseudorep 2
my $full_macs2_pr2_encodePeak_file    = $full_macs2_pr2_basename . '_peaks.encodePeak';
my $full_macs2_pr2_regionPeak_file    = $full_macs2_pr2_basename . '.regionPeak';
my $full_macs2_pr2_regionPeak_gz_file = $full_macs2_pr2_basename . '.regionPeak.gz';

my $full_idr_overlap_file             = $full_idr_basename        . '-overlapped-peaks.txt';
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
my $random_macs2_basename     = $datadir . "\/" . $basename  . '_random_macs2'; #macs2 outputs in these
my $random_macs2_pr1_basename = $datadir . "\/" . $basename  . '_random_macs2_pr1';  #macs2 outputs here for pseudoreplicate 1
my $random_macs2_pr2_basename = $datadir . "\/" . $basename  . '_random_macs2_pr2';  #macs2 outputs here for pseudoreplicate 1
my $random_idr_basename       = $datadir . "\/" . $basename  . '_random_pr1VSpr2';
my $random_idr_plot_basename  = $datadir . "\/" . $basename  . '_random_IDR_plot';#done only for full, not random

#IDR only, MACS2 only
my $random_macs2_encodePeak_file        = $random_macs2_basename     . '_peaks.encodePeak';;
my $random_macs2_regionPeak_file        = $random_macs2_basename     . '.regionPeak';
my $random_macs2_regionPeak_gz_file     = $random_macs2_basename     . '.regionPeak.gz';
#pseudorep 1
my $random_macs2_pr1_encodePeak_file    = $random_macs2_pr1_basename . '_peaks.encodePeak';
my $random_macs2_pr1_regionPeak_file    = $random_macs2_pr1_basename . '.regionPeak';
my $random_macs2_pr1_regionPeak_gz_file = $random_macs2_pr1_basename . '.regionPeak.gz';
#pseudorep 2
my $random_macs2_pr2_encodePeak_file    = $random_macs2_pr2_basename . '_peaks.encodePeak';
my $random_macs2_pr2_regionPeak_file    = $random_macs2_pr2_basename . '.regionPeak';
my $random_macs2_pr2_regionPeak_gz_file = $random_macs2_pr2_basename . '.regionPeak.gz';

my $random_idr_overlap_file             = $random_idr_basename        . '-overlapped-peaks.txt';
my $random_idr_encodepeak_file          = $random_idr_overlap_file    . '.npk'; # generated by idrOverlap2npk.sh
my $temp_idr_encodePeak_file            = $datadir . "\/"             . '_temp.npk';

my $outfile = $datadir . "\/" . $basename . '.result';

#unlink files in case something is left from previous run
unlink $infile_sam;
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
system "$BOWTIE -q -S --best --strata -m 1 -p 20 --chunkmbs 1024 $INDEX_PATH $infile $infile_sam";
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
system "java -Xmx4g -Djava.io.tmpdir=$TMP_DIR -jar $PICARD/SortSam.jar SORT_ORDER=coordinate INPUT=$infile_sam OUTPUT=$infile_bam CREATE_INDEX=true VERBOSITY=ERROR VALIDATION_STRINGENCY=LENIENT";
my $total_mapped_reads = `$SAMTOOLS/samtools view -F 4 $infile_bam | wc -l`; #mapped reads
chomp $total_mapped_reads;
#------------------------------
#bam -> tagalign and bam -> bed (threaded)
#------------------------------
print "-----------------------\n";
print "INFO - create gzipped tagAlign and bed files from initial Bam sample\n";
print "-----------------------\n";
my $thread_bam_to_tagalign  = threads->create( \&bam_to_tagalign, $infile_bam);
my $thread_bam_to_bed       = threads->create( \&bam_to_bed, $infile_bam, $infile_bed);
$thread_bam_to_tagalign->join;
$thread_bam_to_bed->join;
#-------------------------
#peak calling - unfiltered
#-------------------------
print "-----------------------\n";
print "INFO - Peak Calling, full file\n";
print "-----------------------\n";
my $full_thread_gps    = threads->create( \&run_gem,              $infile_bam,         $full_gem_basename    );
my $full_thread_macs14 = threads->create( \&run_macs14,           $infile_tagAlign_gz, $full_macs14_basename );
my $full_thread_macs2  = threads->create( \&run_macs2_get_peaks,  $infile_tagAlign_gz, $full_macs2_basename  );
#genetrack threads, one per chromosome
for my $i (0 .. $#HUMAN_CHRS) {
	push(@full_threads, threads->create(\&run_genetrack, $HUMAN_CHRS[$i], $infile_bed));
}
foreach my $thread (@full_threads){
	my $peak = $thread->join;
	push(@full_peaks, $peak);	
}
my $total_gem_peaks       = $full_thread_gps->join;
my $total_macs14_peaks    = $full_thread_macs14->join;
my $total_macs2_peaks     = $full_thread_macs2->join;
my $total_genetrack_peaks = sum(@full_peaks);
#------------------------
# create pseudoreplicates for MACS2 IDR
#------------------------
print "-----------------------\n";
print "INFO - do_random_split_peakfile.sh: create two random subsets of tagAlign file\n";
print "-----------------------\n";
system "$SCRIPT_PATH/do_random_split_peakfile.sh $infile_tagAlign_gz $datadir";
#-------------------------------
# call peaks on pseudoreplicates (threaded) with MACS2
#-------------------------------
print "-----------------------\n";
print "INFO - MACS2: Call peaks pseudoreplicates\n";
print "-----------------------\n";
my $full_thread_macs2_pr1  = threads->create( \&run_macs2,  $infile_tagAlign_pr1, $full_macs2_pr1_basename);
my $full_thread_macs2_pr2  = threads->create( \&run_macs2,  $infile_tagAlign_pr2, $full_macs2_pr2_basename);
$full_thread_macs2_pr1->join;
$full_thread_macs2_pr2->join;
#---------------------------------
#sort MACS2 output encodePeak file rep1 and rep2 (threaded)
#---------------------------------
print "-----------------------\n";
print "INFO - sort and gzip MACS2 output regionPeak file for pseudorep 1 and 2\n";
print "-----------------------\n";
my $full_thread_sort_pr1  = threads->create( \&sort_pseudorep,  $full_macs2_pr1_encodePeak_file, $full_macs2_pr1_regionPeak_file);
my $full_thread_sort_pr2  = threads->create( \&sort_pseudorep,  $full_macs2_pr2_encodePeak_file, $full_macs2_pr2_regionPeak_file);
$full_thread_sort_pr1->join;
$full_thread_sort_pr2->join;
#------------
#IDR Analysis
#------------
print "-----------------------\n";
print "INFO - IDR run on pseudoreplicates\n";
print "-----------------------\n";
system "$R_CODE/Rscript $IDR/batch-consistency-analysis.r $full_macs2_pr1_regionPeak_file $full_macs2_pr2_regionPeak_file -1 $full_idr_basename 0 F p.value";
#consistency plot? do it only for the full file
system "$R_CODE/Rscript $IDR/batch-consistency-plot.r 1 $full_idr_plot_basename  $full_idr_basename";
my $total_IDR_peaks = `$SCRIPT_PATH/idr_count_peaks.sh  $full_idr_overlap_file`;
chomp $total_IDR_peaks;

#print "MACS2: number of total unfiltered peaks in initial file: "    . $total_macs2_peaks     . "\n";
print "Genetrack: number of total peaks in initial file: " . $total_genetrack_peaks . "\n";
print "GPS: number of total peaks in initial file: "       . $total_gem_peaks       . "\n";
print "MACS14: number of total peaks in initial file: "    . $total_macs14_peaks    . "\n";
print "MACS2: number of total peaks in initial file: "     . $total_macs2_peaks     . "\n";
print "MACS2-IDR: number of total filtered peaks in initial file: "  . $total_IDR_peaks     . "\n";

#remove all temp files
unlink $infile_sam; 
unlink $infile_bai;
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
print $outstream "SAMPLE\tMAPPED_READS\tPEAKS_GENETRACK\tPEAK_RATIO_GENETRACK\tPEAKS_GPS\tPEAK_RATIO_GPS\tPEAKS_MACS14\tPEAK_RATIO_MACS14\tPEAKS_MACS2\tPEAK_RATIO_MACS2\tPEAKS_IDR_MACS2\tPEAK_RATIO_IDR_MACS2\tMEDIAN_ENRICHMENT\n";
system "touch $temp_idr_encodePeak_file";#initally empty, do control first time
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
	print "INFO - create gzipped tagAlign and bed files from random Bam sample\n";
	print "-----------------------\n";
	my $random_thread_bam_to_tagalign  = threads->create( \&bam_to_tagalign, $random_bam_file);
	my $random_thread_bam_to_bed       = threads->create( \&bam_to_bed, $random_bam_file, $random_bed_file);
	$random_thread_bam_to_tagalign->join;
	$random_thread_bam_to_bed->join;
	#-------------------------
	#peak calling - unfiltered
	#-------------------------
	print "-----------------------\n";
	print "INFO - MACS2: Call peaks random tagAlign sample\n";
	print "-----------------------\n";
	my $tmacs14 = threads->create( \&run_macs14,           $random_tagAlign_gz_file, $random_macs14_basename );
	my $tmacs2  = threads->create( \&run_macs2_get_peaks,  $random_tagAlign_gz_file, $random_macs2_basename  );
	my $tgps    = threads->create( \&run_gem,              $random_bam_file, $random_gem_basename            );
	#genetrack threads, one per chromosome
	for my $i (0 .. $#HUMAN_CHRS) {
		push(@gt_threads, threads->create(\&run_genetrack, $HUMAN_CHRS[$i], $random_bed_file));
	}
	#push (@gt_peaks, $_->join) foreach @gt_threads;
	foreach my $thread (@gt_threads){
		my $peak = $thread->join;
		push(@gt_peaks, $peak);	
	}
	my $macs14_peaks = $tmacs14->join;
	my $macs2_peaks  = $tmacs2->join;
	my $gem_peaks    = $tgps->join;
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
	my $random_thread_macs2_pr1  = threads->create( \&run_macs2,  $random_tagAlign_pr1_file, $random_macs2_pr1_basename);
	my $random_thread_macs2_pr2  = threads->create( \&run_macs2,  $random_tagAlign_pr2_file, $random_macs2_pr2_basename);
	$random_thread_macs2_pr1->join;
	$random_thread_macs2_pr2->join;
	#---------------------------------
	#sort MACS2 output encodePeak file rep1 and rep2 (threaded)
	#---------------------------------
	print "-----------------------\n";
	print "INFO - sort and gzip MACS2 output regionPeak file for pseudorep 1 and 2\n";
	print "-----------------------\n";
	my $random_thread_sort_pr1  = threads->create( \&sort_pseudorep,  $random_macs2_pr1_encodePeak_file, $random_macs2_pr1_regionPeak_file);
	my $random_thread_sort_pr2  = threads->create( \&sort_pseudorep,  $random_macs2_pr2_encodePeak_file, $random_macs2_pr2_regionPeak_file);
	$random_thread_sort_pr1->join;
	$random_thread_sort_pr2->join;
	#------------
	#IDR Analysis
	#------------
	print "-----------------------\n";
	print "INFO - IDR run on pseudoreplicates\n";
	print "-----------------------\n";
	system "$R_CODE/Rscript $IDR/batch-consistency-analysis.r $random_macs2_pr1_regionPeak_file $random_macs2_pr2_regionPeak_file -1 $random_idr_basename 0 F p.value";
	my $IDR_peaks = `$SCRIPT_PATH/idr_count_peaks.sh  $random_idr_overlap_file`;
	chomp $IDR_peaks;
	#----------------------
	#convert IDR -overlapped-peaks.txt output in .encodePeak format
	#use slightly modified version of anshul's script
	#USAGE: idrOverlap2npk.sh [idrOverlapFile] [oDir]
	#odir e' opzionale
	#----------------------
	system "$SCRIPT_PATH/idrOverlap2npk.sh $random_idr_overlap_file";
	#output file is $random_idr_encodepeak_file
	#----------------------
	#calculate median fold enrichment using IDR output peaks 
	#----------------------
	#(since there is no control this is the median value in field 7 of the .encodePeak file)
	my $median_newpeak_signal = get_median_peaksignal_for_newpeaks($random_idr_encodepeak_file, $temp_idr_encodePeak_file);
	system "rm $temp_idr_encodePeak_file"; 
	system "cp $random_idr_encodepeak_file  $temp_idr_encodePeak_file"; #copy encode peak in temp file for next iteration
	
	print "MACS2: number of unfiltered peaks in sample "       . $counter . ": "  . $macs2_peaks           . "\n";
	print "IDR: number of total IDR filtered peaks in sample " . $counter . ": "  . $IDR_peaks             . "\n";
	print "Median New peak signal "                            . $counter . ": "  . $median_newpeak_signal . "\n";
	
	print $outstream  $counter 
	                  . "\t" . $number_of_mapped_reads
	                  . "\t" . $genetrack_peaks  . "\t" . ($genetrack_peaks/$total_genetrack_peaks)  
	                  . "\t" . $gem_peaks        . "\t" . ($gem_peaks/$total_gem_peaks) 
	                  . "\t" . $macs14_peaks     . "\t" . ($macs14_peaks/$total_macs14_peaks)
					  . "\t" . $macs2_peaks      . "\t" . ($macs2_peaks/$total_macs2_peaks) 
					  . "\t" . $IDR_peaks        . "\t" . ($IDR_peaks/$total_IDR_peaks)
					  . "\t" . $median_newpeak_signal
					  . "\n" ;

	$BIN += $INITIAL_BIN;
	$counter += 1;
	unlink $random_bam_file;
	unlink $random_bai_file;
	unlink $random_tagAlign_gz_file;
	unlink $random_tagAlign_pr1_file;
	unlink $random_tagAlign_pr2_file;

	system "rm $gt_outputs_gff"; #rm genetrack outputs
	system "rm $gt_outputs_bed"; #rm genetrack outputs
    system "rm $random_macs14_basename*.*"; #rm macs14 output
	system "rm $random_macs2_basename*.*"; #rm macs2 output
	system "rm -rf $random_idr_basename*"; #remove IDR output
}
print $outstream  $counter 
                  . "\t" . $total_mapped_reads
                  . "\t" . $total_genetrack_peaks  . "\t" . ($total_genetrack_peaks/$total_genetrack_peaks)  
                  . "\t" . $total_gem_peaks        . "\t" . ($total_gem_peaks/$total_gem_peaks) 
                  . "\t" . $total_macs14_peaks     . "\t" . ($total_macs14_peaks/$total_macs14_peaks)
				  . "\t" . $total_macs2_peaks      . "\t" . ($total_macs2_peaks/$total_macs2_peaks) 
				  . "\t" . $total_IDR_peaks        . "\t" . ($total_IDR_peaks/$total_IDR_peaks)
				  . "\t" . ''
				  . "\n" ;
unlink $infile_bam;
close $outstream;
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
	my  ( $encodePeak_file, $previous_encodePeak_file ) = @_;
	my $working_peak_file = $datadir . "\/" .  '_working.npk'; ;
	system "rm $working_peak_file";
	my @novel_peak_signals;

	if (-s $previous_encodePeak_file) { #if the temp file is not empty
		system "$BEDTOOLS/subtractBed -A -a $encodePeak_file -b $previous_encodePeak_file  > $working_peak_file";
	}else{
		system "cp $encodePeak_file $working_peak_file";#it means this is the first iteration	
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

	$stat->add_data(@novel_peak_signals);
	print $stat->mean(), " is the mean", "\n";
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
	my ($bam_input, $bed_output);
	system "$BEDTOOLS/bamToBed -i $bam_input > $bed_output";
	if ( $? == -1 ){
    	print "bam_to_bed(): problem with output: $!\n";
    	exit -1;
	}
	return;
}
#-----------------------------------------------------------------------------
sub sort_pseudorep{
	my ( $encodePeak_file, $regionPeak_file ) = @_;
	system "sort -k 8nr,8nr $encodePeak_file | head -n 100000 > $regionPeak_file";
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
	my $temp_file              = $datadir . "\/" . $input_id .  '_'  . $chr . '_' .  'random_raw.gff';
	my $random_gff_file        = $datadir . "\/" . $input_id .  '_'  . $chr . '_' .  'random.gff';
	my $random_sorted_gff_file = $datadir . "\/" . $input_id .  '_'  . $chr . '_' .  'random_sorted.gff';
	my $random_sorted_bed_file = $datadir . "\/" . $input_id .  '_'  . $chr . '_' .  'random.bed'; # this is the one where you need to count

	#python ${PCODE}/genetrack.py -v -s 5 -e 10 -F 6 -c chr1 ${FILE}
	system "python $GENETRACK -s 5 -e 10 -c $chr $input_file > $temp_file";
    if ( $? == -1 ){
    	print "Genetrack: problem with output: $!\n";
        exit -1;
    }
	#1 - remove all "singletons" (entries where sd = 0.0)
	open (my $temp_data,  q{<}, $temp_file) or die ("Unable to open $temp_file: $!");
	open (my $out_data,  q{>}, $random_gff_file) or die ("Unable to open file $random_gff_file: $!");
	while(<$temp_data>){
		print $out_data $_ unless($_ =~ /stddev=0.0;/);
	}
	close $temp_data;
	close $out_data;
	#2 - merge + and - peaks
	#sort by chr and by start coordinate first (required by bedtools)
	system "sort -k1,1 -k4n,4 $random_gff_file > $random_sorted_gff_file";
	if ( $? == -1 ){
		print "sort: problem with output: $!\n";
		exit -1;
	}
	#YOU NEED TO OPTIMISE THAT -d!!
	system "$BEDTOOLS/mergeBed -n -d 35 -i $random_sorted_gff_file > $random_sorted_bed_file";
	if ( $? == -1 ){
		print "mergeBed: problem with output: $!\n";
		exit -1;
	}
	my $wc_out = `cat $random_sorted_bed_file | wc -l`;
	return (chomp $wc_out);
}
#-----------------------------------------------------------------------------
sub run_gem{
	my ( $input_file, $gem_basename ) = @_;
	#system "java -Xmx10G -jar $GEM/gem.jar --g $CHROM_SIZE_PATH --d $GEM/Read_Distribution_ChIP-exo.txt --s $MAPPABILITY --expt $random_bam_file --f SAM --outBED --out $random_gem_basename --genome $GENOME_PATH --k_min 2 --k_max 40";
	system "java -Xmx10G -jar $GEM/gem.jar --g $CHROM_SIZE_PATH --d $GEM/Read_Distribution_ChIP-exo.txt --s $MAPPABILITY --expt $input_file --f SAM --top 0 --out $gem_basename";
	if ( $? == -1 ){
                print "GEM: problem with output: $!\n";
                exit -1;
    }
    #now capture number of rows in bedfile and fill an output file with them
    my $gem_peakfile = $gem_basename . "_GPS_events.txt";
    my $wc_out = `cat $gem_peakfile | wc -l`;
	return (chomp $wc_out);
}
#-----------------------------------------------------------------------------
sub run_macs14{
	my ( $input_file, $macs14_basename ) = @_;
	#/macs14 -t ${FILE} -n ${PDATA}/${ID} -s 40 --nomodel --shiftsize=13 --keep-dup=all -p 1e-8 -g 2540757438
	#system "macs14 -t  $input_file -n $macs14_basename -g $MAPPABILITY -s 40 --keep-dup=all --nomodel --shiftsize 13 --verbose 1 --nolambda --llocal 0";
	system "macs14 -t  $input_file -f BED -n $macs14_basename -g $MAPPABILITY -s 40 --keep-dup=all --nomodel --shiftsize 13 --verbose 1";	
	if ( $? == -1 ){
		print "MACS14: problem with output: $!\n";
                exit -1;
        }
	#now capture number of rows in bedfile and fill an outputfile with the
	my $macs14_bedfile = $macs14_basename . "_peaks.bed";
    my $wc_out = `cat $macs14_bedfile | wc -l`;
	return (chomp $wc_out);	
}
#-----------------------------------------------------------------------------
sub run_macs2_get_peaks{
	my ( $input_file, $macs2_basename ) = @_;
	
	run_macs2($input_file, $macs2_basename);
	#now capture number of rows in narrowPeak bed file
	my $macs2_bedfile = $macs2_basename . "_peaks.encodePeak";
    my $wc_out = `cat $macs2_bedfile | wc -l`;
	return (chomp $wc_out);
}
#-----------------------------------------------------------------------------
sub run_macs2{
	my ( $input_file, $macs2_basename ) = @_;
	#${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/${ID} -g 2540757438 --keep-dup=all -q 0.001 --nomodel --extsize 26 --call-summits
	#system "macs2 callpeak -t $input_file -n $macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 --verbose 1 --nolambda --llocal 0";
	#use the following from IDR:
	#macs2 callpeak -t $infile_tagAlign_gz -f BED -n $full_macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 -p 1e-3 --to-large
	#system "macs2 callpeak -t $input_file -f BED -n $macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 -p 1e-3 --verbose 1 --nolambda --llocal 0";
	system "macs2 callpeak -t $input_file -f BED -n $macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 --verbose 1";
	if ( $? == -1 ){
		print "MACS2: problem with output: $!\n";
        exit -1;
    }
    return;	
}