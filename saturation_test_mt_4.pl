#!/usr/bin/perl
use strict;
use warnings;
use threads;
use Getopt::Long;
use File::Basename;
use List::Util qw(sum);

#13/3/2013
#This is obtained from saturation_test_mt_2.pl
#With the following differences
#only do macs2
#do the full postprocessing IDR pipeline (standalone in the _idr_pseudorep.pl file) after each MACS2 peak calling
#only count the peaks after IDR thresholding

#IDR is very VERY slow, hence why doing only MACS2
#also Anshul suggests not to use on MACS14
#should probably also try it with GEM
#should probably parallelise the random subsets, for now they happen in sequence and all output for one iteration is deleted after output peak number is collecred

#INPUTS
#1 fastq file
#2 step or frequency (0.01? 0.1)

#OUTPUTS
#TSV FILE  <uniquely mapped reads> <number of peaks> <number of peaks/total peaks for full file>

my $infile;
my $STEP;
my $input_id; my $datadir; # to customise randomised outputs
my $SCRIPT_PATH = '/net/isi-scratch/giuseppe/scripts';
my $TOOL_PATH = '/net/isi-scratch/giuseppe/tools';
my $INDEX_PATH = '/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19_bowtie';
my $CHROM_SIZE_PATH = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes';
my $GENOME_PATH = '/net/isi-scratch/giuseppe/indexes/Hsap/hg19_genome';
my $R_CODE = '/net/isi-cgat/ifs/apps/apps/R-2.14.1/bin';
my $TMP_DIR = '/tmp';

#program paths, change as needed
my $BOWTIE   = $TOOL_PATH . '/bowtie-0.12.9/bowtie';
my $PICARD   = $TOOL_PATH . '/picard-tools-1.84/picard-tools-1.84';
my $SAMTOOLS = $TOOL_PATH . '/samtools';
my $BEDTOOLS = $TOOL_PATH . '/bedtools-2.17.0/bin';
my $IDR = $TOOL_PATH . '/idrCode'; #now using idr in giuseppe/scripts
#my $SPP = $TOOL_PATH . '/phantompeakqualtools';

#human- and dataset-specific
my $MAPPABILITY = 2540757438; #for Hsap (40bp) used by Macs2 and GEM

GetOptions(
	'input=s'  =>\$infile,
	'step=f' => \$STEP,
);
#constants
my $INIT = $STEP; #maybe change this to something different
my $MAX = 1;

if(!$infile){
     print "USAGE: do_saturation_test.pl -input=<INFILE> -step=<STEP>\n";
     print "<INFILE> input fastq file\n";
     print "<STEP>: increment for additive randomiser fraction of reads (e.g. 0.01)\n";
     exit 1;
}

if(!$STEP){
     print "USAGE: do_saturation_test.pl -input=<INFILE> -step=<STEP>\n";
     print "<INFILE> input fastq file\n";
     print "<STEP>: increment for additive randomiser fraction of reads (e.g. 0.01)\n";
     exit 1;
}


#get unique id to customise the random outputs
#TODO this is heavily dependent on the input files and will work only with VDR inputs now
if($infile =~ /(VDR\_)(\w{2}\d{5})(\_Peconic\d{5})(\_trimmed)(.fastq)$/){  #eg VDR_GM19213_Peconic20324_trimmed.fastq
	$input_id = $2;
}else{
	print "Input file name not recognised. Modify the regular expression in the script.\n";
	exit -1;
}
#I want all the outputs in a directory whose format will be $PATH/$input_id_saturation_data
my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
if($directory eq "\.\/"){ 
	$datadir = $input_id . '_saturation_test';
}else{
	$datadir =  $directory . "\/" . $input_id . '_saturation_test';
}
system "mkdir $datadir";

#---------
#data - complete
#---------
my $infile_sam              = $datadir . "\/" . $basename . '.sam';
my $infile_bam              = $datadir . "\/" . $basename . '.bam';#this is now kept throughout the whole duration of the script; picard will sample from this
my $infile_bai              = $datadir . "\/" . $basename . '.bai';
my $infile_tagAlign         = $datadir . "\/" . $basename . '.tagAlign';
my $infile_tagAlign_gz      = $datadir . "\/" . $basename . '.tagAlign.gz';
my $infile_tagAlign_pr1     = $datadir . "\/" . $basename . '.pr1.tagAlign.gz';
my $infile_tagAlign_pr2     = $datadir . "\/" . $basename . '.pr2.tagAlign.gz';

my $full_macs2_basename     = $datadir . "\/" . $basename  . '_macs2';      #macs2 outputs in these
my $full_macs2_pr1_basename = $datadir . "\/" . $basename  . '_macs2_pr1';  #macs2 outputs here for pseudoreplicate 1
my $full_macs2_pr2_basename = $datadir . "\/" . $basename  . '_macs2_pr2';  #macs2 outputs here for pseudoreplicate 1
my $full_idr_basename       = $datadir . "\/" . $basename  . '_pr1VSpr2';
my $full_idr_plot_basename  = $datadir . "\/" . $basename  . '_IDR_plot';#done only for full, not random

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
my $random_tagAlign_file      = $datadir . "\/" . $basename  . '_random.tagAlign';
my $random_tagAlign_gz_file   = $datadir . "\/" . $basename  . '_random.tagAlign.gz';
my $random_tagAlign_pr1_file  = $datadir . "\/" . $basename  . '_random.pr1.tagAlign.gz';
my $random_tagAlign_pr2_file  = $datadir . "\/" . $basename  . '_random.pr2.tagAlign.gz';

my $random_macs2_basename     = $datadir . "\/" . $basename  . '_random_macs2'; #macs2 outputs in these
my $random_macs2_pr1_basename = $datadir . "\/" . $basename  . '_random_macs2_pr1';  #macs2 outputs here for pseudoreplicate 1
my $random_macs2_pr2_basename = $datadir . "\/" . $basename  . '_random_macs2_pr2';  #macs2 outputs here for pseudoreplicate 1
my $random_idr_basename       = $datadir . "\/" . $basename  . '_random_pr1VSpr2';
my $random_idr_plot_basename  = $datadir . "\/" . $basename  . '_random_IDR_plot';#done only for full, not random

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

my $outfile = $datadir . "\/" . $basename . '.result';

#unlink files in case something is left from previous run
unlink $infile_sam;
unlink $infile_bam;
unlink $infile_bai;
unlink $infile_tagAlign_gz;
unlink $infile_tagAlign_pr1;
unlink $infile_tagAlign_pr2;
unlink $random_bam_file;
unlink $random_bai_file;
unlink $random_tagAlign_gz_file;
unlink $random_tagAlign_pr1_file;
unlink $random_tagAlign_pr2_file;

#############################################################
#run pipeline on full dataset to get total number of peaks
#############################################################
my $total_reads = `wc -l $infile`;
my @total_reads = split(/ /, $total_reads);
my $total_fq_reads = $total_reads[0]/4;	
#------------------
#align (only once)
#------------------
print "Aligning initial .fastq file..\n";
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
my $total_mapped_reads = `$SAMTOOLS/samtools view -F 4 $infile_bam | wc -l`; #uniquely mapped reads
chomp $total_mapped_reads;

#----------------
#bam -> tagalign
#----------------
print "-----------------------\n";
print "INFO - do_bam_to_tagalign.sh: create gzipped tagAlign from initial Bam file\n";
print "-----------------------\n";
system "$SCRIPT_PATH/do_bam_to_tagalign.sh $infile_bam $datadir";
#-------------------------
#peak calling - unfiltered
#-------------------------
print "-----------------------\n";
print "INFO - MACS2: Call peaks full file\n";
print "-----------------------\n";
my $total_macs2_peaks = run_macs2_get_peaks($infile_tagAlign_gz, $full_macs2_basename);
#------------------------
# create pseudoreplicates
#------------------------
print "-----------------------\n";
print "INFO - do_random_split_peakfile.sh: create two random subsets of tagAlign file\n";
print "-----------------------\n";
system "$SCRIPT_PATH/do_random_split_peakfile.sh $infile_tagAlign_gz $datadir";

#-------------------------------
# call peaks on pseudoreplicates (threaded)
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

print "MACS2: number of total unfiltered peaks in initial file: "    . $total_macs2_peaks     . "\n";
print "MACS2: number of total IDR filtered peaks in initial file: "  . $total_IDR_peaks     . "\n";


exit;

#remove all temp files
unlink $infile_sam; 
unlink $infile_bai;
unlink $infile_tagAlign_gz;
unlink $infile_tagAlign_pr1;
unlink $infile_tagAlign_pr2;
system "rm -rf $full_macs2_basename*"; #remove mac2 output
system "rm -rf $full_idr_basename*"; #remove IDR output (but not the plot)

#########################
#randomisation main loop
#########################
#outstream header
open (my $outstream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
print $outstream "SAMPLE\tMAPPED_READS_BAM\tPEAKS_MACS2\tPEAK_RATIO_MACS2\tPEAKS_IDR_MACS2\tPEAK_RATIO_IDR_MACS2\n";
my $counter = 1;
while ($INIT <= $MAX){
	my $wc_out; my @wc_out;
	#----------------
	#random subset selection from bam file (Picard)
	#----------------
	print "Creating random .bam file, (fraction: $INIT)..\n";
    system "java -Xmx4g -Djava.io.tmpdir=$TMP_DIR -jar $PICARD/DownsampleSam.jar INPUT=$infile_bam OUTPUT=$random_bam_file RANDOM_SEED=null VERBOSITY=ERROR PROBABILITY=$INIT";
    my $number_of_mapped_reads = `$SAMTOOLS/samtools view -F 4 $random_bam_file | wc -l`; #uniquely mapped reads
	chomp $number_of_mapped_reads;
	#----------------
	#bam -> tagalign
	#----------------
	print "-----------------------\n";
	print "INFO - do_bam_to_tagalign.sh: create gzipped tagAlign from random Bam sample\n";
	print "-----------------------\n";
	system "$SCRIPT_PATH/do_bam_to_tagalign.sh $random_bam_file $datadir";
	#-------------------------
	#peak calling - unfiltered
	#-------------------------
	print "-----------------------\n";
	print "INFO - MACS2: Call peaks random tagAlign sample\n";
	print "-----------------------\n";
	my $random_macs2_peaks = run_macs2_get_peaks($random_tagAlign_gz_file, $random_macs2_basename);
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
	my $random_IDR_peaks = `$SCRIPT_PATH/idr_count_peaks.sh  $random_idr_overlap_file`;
	chomp $random_IDR_peaks;
	
	print "MACS2: number of unfiltered peaks in sample "       . $counter . ": "  . $random_macs2_peaks     . "\n";
	print "IDR: number of total IDR filtered peaks in sample " . $counter . ": "  . $random_IDR_peaks  . "\n";
	
	print $outstream  $counter . "\t" .  $number_of_mapped_reads . "\t" . $random_macs2_peaks . "\t" . ($random_macs2_peaks/$total_macs2_peaks) . "\t" . $random_IDR_peaks . "\t" . ($random_IDR_peaks/$total_IDR_peaks) . "\n" ;

	$INIT += $STEP;
	$counter += 1;
	unlink $random_bam_file;
	unlink $random_bai_file;
	unlink $random_tagAlign_gz_file;
	unlink $random_tagAlign_pr1_file;
	unlink $random_tagAlign_pr2_file;

	system "rm $random_macs2_basename*.*"; #rm macs2 output
	system "rm -rf $random_idr_basename*"; #remove IDR output
}
print $outstream  $counter . "\t" . $total_mapped_reads . "\t" .  $total_macs2_peaks . "\t" . ($total_macs2_peaks/$total_macs2_peaks) . "\t" . $total_IDR_peaks . "\t" . ($total_IDR_peaks/$total_IDR_peaks) .  "\n";
unlink $infile_bam;
close $outstream;
print "\nFINISHED\n";


sub sort_pseudorep{
	my ( $encodePeak_file, $regionPeak_file ) = @_;
	system "sort -k 8nr,8nr $encodePeak_file | head -n 100000 > $regionPeak_file";
	return;
}

sub run_macs2_get_peaks{
	my ( $input_file, $macs2_basename ) = @_;
	
	run_macs2($input_file, $macs2_basename);
	#now capture number of rows in bedfile and fill an outputfile with the
	my $macs2_bedfile = $macs2_basename . "_peaks.encodePeak";
    my $wc_out = `wc -l $macs2_bedfile`;
    my @peaks = split(/ /, $wc_out);
	return $peaks[0];	
}


#basic peak calling routine for parallelisation in pseudoreplicate count
#  --nolambda            If True, MACS will use fixed background lambda as
#                        local lambda for every peak region. Normally, MACS
#                        calculates a dynamic local lambda to reflect the local
#                        bias due to potential chromatin structure.

# --llocal LARGELOCAL   The large nearby region in basepairs to calculate
#                        dynamic lambda. This is used to capture the surround
#                        bias. If you set this to 0, MACS will skip llocal
#                        lambda calculation. *Note* that MACS will always
#                        perform a d-size local lambda calculation. The final
#                        local bias should be the maximum of the lambda value
#                        from d, slocal, and llocal size windows. DEFAULT:
#                        10000.
sub run_macs2{
	my ( $input_file, $macs2_basename ) = @_;
	#${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/${ID} -g 2540757438 --keep-dup=all -q 0.001 --nomodel --extsize 26 --call-summits
	#system "macs2 callpeak -t $input_file -n $macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 --verbose 1 --nolambda --llocal 0";
	#use the following from IDR:
	#macs2 callpeak -t $infile_tagAlign_gz -f BED -n $full_macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 -p 1e-3 --to-large
	system "macs2 callpeak -t $input_file -f BED -n $macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 -p 1e-3 --verbose 1 --nolambda --llocal 0";
	if ( $? == -1 ){
		print "MACS2: problem with output: $!\n";
        exit -1;
    }
    return;	
}






