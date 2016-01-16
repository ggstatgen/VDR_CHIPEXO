#!/usr/bin/perl
use strict;
use warnings;
use threads;
use Getopt::Long;
use File::Basename;
use List::Util qw(sum);

#7/3/2013
#Inspired from the Landt 2012 paper.
#Calculate the fraction of all mapped reads that fall into peak regions
#PROCEDURE
#Read peak file from peak caller output (bed? xls? genetrack before bedtools conversion?)
##need number of uniquely mapped reads
#count number of uniquely mapped reads falling under peaks

#do FRiP values correlate positively and linearly with number of peaks?
#as you call more peaks, do you see higher FRiP?


#INPUTS
#1 fastq file
#2 step or frequency (0.01? 0.1)

#OUTPUTS
#TSV FILE  <uniquely mapped reads> <number of peaks> <number of peaks/total peaks for full file>
my $infile;
my $input_id; my $datadir; # to customise randomised outputs
my $SCRIPT_PATH = '/net/isi-scratch/giuseppe/scripts';
my $CODE_PATH = '/net/isi-scratch/giuseppe/tools';
my $INDEX_PATH = '/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19_bowtie';
my $CHROM_SIZE_PATH = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes';
my $GENOME_PATH = '/net/isi-scratch/giuseppe/indexes/Hsap/hg19_genome';
my $TMP_DIR = '/tmp';

#program paths, change as needed
my $BOWTIE = $CODE_PATH . '/bowtie-0.12.9/bowtie';
my $PICARD = $CODE_PATH . '/picard-tools-1.84/picard-tools-1.84';
my $GENETRACK = $CODE_PATH . '/chipexo-master/genetrack/genetrack.py';
my $SAMTOOLS = $CODE_PATH . '/samtools';
my $BEDTOOLS = $CODE_PATH . '/bedtools-2.17.0/bin';
my $GEM = $CODE_PATH . '/gem';

my $STEP;
#constants
my $INIT = 0.01; #maybe change this to something different
my $MAX = 1;

#human- and dataset-specific
my $MAPPABILITY = 2540757438; #for Hsap (40bp) used by Macs2 and GEM
my @HUMAN_CHRS = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'); 

GetOptions(
	'input=s'  =>\$infile,
	'step=f' => \$STEP,
);

if(!$infile){
     print "USAGE: FRiP_test.pl -input=<INFILE> -step=<STEP>\n";
     print "<INFILE> input fastq file\n";
     print "<STEP>: increment for additive randomiser fraction of reads (e.g. 0.01)\n";
     exit 1;
}

if(!$STEP){
     print "USAGE: FRiP_test.pl -input=<INFILE> -step=<STEP>\n";
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
#I want all the outputs in a directory whose format will be $PATH/$input_id_FRiP_data
my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
if($directory eq "\.\/"){ 
	$datadir = $input_id . '_FRiP_test';
}else{
	$datadir =  $directory . "\/" . $input_id . '_FRiP_test';
}
system "mkdir $datadir";
#---------
#temp data
#---------
my $infile_sam = $datadir . "\/" . $basename . '.sam';
my $infile_bam = $datadir . "\/" . $basename . '.bam';#this is now kept throughout the whole duration of the script; picard will sample from this
my $infile_bai = $datadir . "\/" . $basename . '.bai';
my $infile_bed = $datadir . "\/" . $basename . '.bed';
my $full_gem_basename        = $datadir . "\/" . $basename  . '_gem'; #gem outputs in these
my $full_macs14_basename     = $datadir . "\/" . $basename  . '_macs14'; #macs14 outputs in these
my $full_macs2_basename      = $datadir . "\/" . $basename  . '_macs2';  #macs2 outputs in these

my $outfile = $datadir . "\/" . $basename . '.result';

my $random_sam_file         = $datadir . "\/" . $input_id  . '_random.sam'; #bowtie writes here / picard reads
my $random_bam_file         = $datadir . "\/" . $input_id  . '_random.bam'; #picard writes here /macs2 and GEM read
my $random_bai_file         = $datadir . "\/" . $input_id  . '_random.bai'; #picard writes here
my $random_bed_file         = $datadir . "\/" . $input_id  . '_random.bed'; #bamtobed writes here / genetrack reads
my $random_macs2_basename   = $datadir . "\/" . $input_id  . '_random_macs2'; #macs2 outputs in these
my $random_macs14_basename  = $datadir . "\/" . $input_id  . '_random_macs14'; #macs14 outputs in these
my $random_gem_basename     = $datadir . "\/" . $input_id  . '_random_gem'; #gem outputs in these

#unlink files in case something is left from previous run
unlink $random_sam_file; 
unlink $random_bam_file;
unlink $random_bai_file;
unlink $random_bed_file;
unlink $infile_sam;
unlink $infile_bam;
unlink $infile_bed;

#############################################################
#run pipeline on full dataset to get total number of peaks
#############################################################
my @full_threads; my @full_peaks;
my $total_reads = `wc -l $infile`;
my @total_reads = split(/ /, $total_reads);
my $total_fq_reads = $total_reads[0]/4;	
#--------------
#align
#-------------
print "Aligning initial .fastq file..\n";
system "$BOWTIE -q -S --best --strata -m 1 -p 12 --chunkmbs 1024 $INDEX_PATH $infile $infile_sam";
if ( $? == -1 ){
	print "BOWTIE: problem with output: $!\n";
	exit -1;
}
#-------------
#sam to bam
#------------
system "java -Xmx4g -Djava.io.tmpdir=$TMP_DIR -jar $PICARD/SortSam.jar SORT_ORDER=coordinate INPUT=$infile_sam OUTPUT=$infile_bam CREATE_INDEX=true VERBOSITY=WARNING VALIDATION_STRINGENCY=LENIENT";
my $total_mapped_reads = `$SAMTOOLS/samtools view -F 4 $infile_bam | wc -l`; #uniquely mapped reads
chomp $total_mapped_reads;
#------------
#bam to bed (only needed by genetrack threads - improve multithreading here)
#------------
#convert bam to bed
#${PCODE}/bamToBed -i ${FILE} > ${PDATA}/${BASEFILE}.bed;
#system "$BEDTOOLS/bamToBed -i $infile_bam > $infile_bed";
#if ( $? == -1 ){
#	print "bamToBed: problem with output: $!\n";
#    exit -1;
#}
#---------------
#peak calling
#---------------
#gem peaks
my $full_thread_gps    = threads->create( \&run_gem,    $infile_bam, $full_gem_basename );
my $full_thread_macs14 = threads->create( \&run_macs14, $infile_bam, $full_macs14_basename );
my $full_thread_macs2  = threads->create( \&run_macs2,  $infile_bam, $full_macs2_basename  );
#genetrack threads, one per chromosome
#for my $i (0 .. $#HUMAN_CHRS) {
#	push(@full_threads, threads->create(\&run_genetrack, $HUMAN_CHRS[$i], $infile_bed));
#}
#foreach my $thread (@full_threads){
#	my $peak = $thread->join;
#	push(@full_peaks, $peak);	
#}
my $total_gem_peaks    = $full_thread_gps->join;
my $total_macs14_peaks = $full_thread_macs14->join;
my $total_macs2_peaks  = $full_thread_macs2->join;
#my $total_genetrack_peaks = sum(@full_peaks);

#print "Genetrack: number of total peaks in initial file: " . $total_genetrack_peaks . "\n";
print "GPS: number of total peaks in initial file: "       . $total_gem_peaks       . "\n";
print "MACS14: number of total peaks in initial file: "    . $total_macs14_peaks    . "\n";
print "MACS2: number of total peaks in initial file: "     . $total_macs2_peaks     . "\n";

#---------------
#compute FRiP: input: total unique mapped reads, xsl output of macs
#---------------
my $total_FRiP_macs14 = do_FRiP_macs($total_mapped_reads, $full_macs14_basename);
my $total_FRiP_macs2 = do_FRiP_macs($total_mapped_reads, $full_macs2_basename);
my $total_FRiP_gem = do_FRiP_gem($total_mapped_reads, $full_gem_basename);

print $total_mapped_reads, "\n"; 
print 'FRIPS 1.4: ', $total_FRiP_macs14, "\n";
print 'FRIPS 2: ', $total_FRiP_macs2, "\n";
print 'FRIPS gem: ', $total_FRiP_gem, "\n";

#remove all temp files
unlink $infile_sam;
#unlink $infile_sam;#this is needed in the randomisation lopp
unlink $infile_bai;
unlink $infile_bed;
#my $gt_outputs_gff = $datadir . "\/" . "$input_id\*chr\*.gff";
#my $gt_outputs_bed = $datadir . "\/" . "$input_id\*chr\*.bed";
#system "rm $gt_outputs_gff"; #rm genetrack outputs
#system "rm $gt_outputs_bed"; #rm genetrack outputs
system "rm -rf $full_gem_basename*"; #remove gem output
system "rm -rf $full_macs14_basename*"; #remove macs14 output
system "rm -rf $full_macs2_basename*"; #remove macs2 output



#########################
#randomisation main loop
#########################
#outstream header
open (my $outstream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
print $outstream "SAMPLE\tMAPPED_READS_BAM\tPEAKS_GEM\tFRIP_GEM\tPEAKS_MACS14\tFRIP_MACS14\tPEAKS_MACS2\tFRIPS_MACS2\n";
my $counter = 1;
while ($INIT <= $MAX){
	my $wc_out; my @wc_out;
	my @gt_threads; my @gt_peaks;
	#----------------
	#random subset selection from bam file (Picard)
	#----------------
	print "Creating random .bam file, (fraction: $INIT)..\n";
    system "java -Xmx4g -Djava.io.tmpdir=$TMP_DIR -jar $PICARD/DownsampleSam.jar INPUT=$infile_bam OUTPUT=$random_bam_file RANDOM_SEED=null VERBOSITY=INFO PROBABILITY=$INIT";
    my $number_of_mapped_reads = `$SAMTOOLS/samtools view -F 4 $random_bam_file | wc -l`; #uniquely mapped reads
	chomp $number_of_mapped_reads;
	#------------
	#bam to bed (only needed by genetrack threads - improve )
	#------------
	#${PCODE}/bamToBed -i ${FILE} > ${PDATA}/${BASEFILE}.bed;
#	system "$BEDTOOLS/bamToBed -i $random_bam_file > $random_bed_file";
#	if ( $? == -1 ){
#                print "bamToBed: problem with output: $!\n";
#                exit -1;
#	}
	#---------------
	#peak calling
	#---------------
	my $tmacs14 = threads->create( \&run_macs14, $random_bam_file, $random_macs14_basename );
	my $tmacs2  = threads->create( \&run_macs2,  $random_bam_file, $random_macs2_basename  );
	my $tgps    = threads->create( \&run_gem,    $random_bam_file, $random_gem_basename    );
	#genetrack threads, one per chromosome
#	for my $i (0 .. $#HUMAN_CHRS) {
#		push(@gt_threads, threads->create(\&run_genetrack, $HUMAN_CHRS[$i], $random_bed_file));
#	}
#	#push (@gt_peaks, $_->join) foreach @gt_threads;
#	foreach my $thread (@gt_threads){
#		my $peak = $thread->join;
#		push(@gt_peaks, $peak);	
#	}
	my $macs14_peaks = $tmacs14->join;
	my $macs2_peaks  = $tmacs2->join;
	#my $genetrack_peaks = sum(@gt_peaks);
	my $gem_peaks = $tgps->join;
	
	my $random_FRiP_macs14 = do_FRiP_macs($number_of_mapped_reads, $random_macs14_basename);
	my $random_FRiP_macs2  = do_FRiP_macs($number_of_mapped_reads, $random_macs2_basename);
	my $random_FRiP_gem    = do_FRiP_gem($number_of_mapped_reads, $random_gem_basename);

	print $outstream  $counter . "\t" .  $number_of_mapped_reads . "\t" . $gem_peaks . "\t" . $random_FRiP_gem . "\t" . $macs14_peaks . "\t" . $random_FRiP_macs14 .  "\t" . $macs2_peaks . "\t" . $random_FRiP_macs2 . "\n" ;

	$INIT += $STEP;
	$counter += 1;
	unlink $random_sam_file;
	unlink $random_bam_file;
	unlink $random_bai_file;
	unlink $random_bed_file;

   	#unlink $random_gff_file_s;
	#system "rm $gt_outputs_gff"; #rm genetrack outputs
	#system "rm $gt_outputs_bed"; #rm genetrack outputs
	system "rm $random_macs2_basename*.*"; #rm macs2 output
	system "rm $random_macs14_basename*.*"; #rm macs14 output
	#system "rm -rf $random_gem_basename*"; #remove gem output
}
print $outstream  $counter . "\t" . $total_mapped_reads . "\t" . $total_gem_peaks . "\t" . $total_FRiP_gem . "\t" . $total_macs14_peaks . "\t" . $total_FRiP_macs14  .  "\t" . $total_macs2_peaks . "\t" . $total_FRiP_macs2 . "\n";
unlink $infile_bam;
close $outstream;


sub do_FRiP_gem{
	my ( $total_reads, $basename ) = @_;
	#open even file
	my $peak_data = $basename . "_GPS_events.txt";
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
	print 'GEM: Total reads in peaks: ', $total_reads_in_peaks, "\n";
	return ($total_reads_in_peaks/$total_reads);
}


#this gets an .xls file produced my macs, and the total number or mapped reads
#for eack peak in the xls, it gets the "pileup" column, sums for all peaks
#it returns sum(pileup)/(total number of aligned reads)
sub do_FRiP_macs{
	my ( $total_reads, $basename ) = @_;
	#open excel file
	my $peak_data = $basename . '_peaks.xls'; #macs14 and macs2 only
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
	print 'MACS: Total reads in peaks: ', $total_reads_in_peaks, "\n";
	return ($total_reads_in_peaks/$total_reads);
}



#genetrack output is as follows
#chr1    genetrack       .       568373  568383  6.24515535516   +       .       readcount=6;ID=P568378;stddev=2.35702260396;height=6.24515535516
#chr1    genetrack       .       568477  568487  6.0292452985    +       .       readcount=5;ID=P568482;stddev=0.0;height=6.0292452985
#sub run_genetrack{
#	my ( $chr, $input_file ) = @_; 
#	my $temp_file              = $datadir . "\/" . $input_id .  '_'  . $chr . '_' .  'random_raw.gff';
#	my $random_gff_file        = $datadir . "\/" . $input_id .  '_'  . $chr . '_' .  'random.gff';
#	my $random_sorted_gff_file = $datadir . "\/" . $input_id .  '_'  . $chr . '_' .  'random_sorted.gff';
#	my $random_sorted_bed_file = $datadir . "\/" . $input_id .  '_'  . $chr . '_' .  'random.bed'; # this is the one where you need to count
#
#	#python ${PCODE}/genetrack.py -v -s 5 -e 10 -F 6 -c chr1 ${FILE}
#	system "python $GENETRACK -s 5 -e 10 -c $chr $input_file > $temp_file";
#    if ( $? == -1 ){
#    	print "Genetrack: problem with output: $!\n";
#        exit -1;
#    }
#	
#	#1 - remove all "singletons" (entries where sd = 0.0)
#	open (my $temp_data,  q{<}, $temp_file) or die ("Unable to open $temp_file: $!");
#	open (my $out_data,  q{>}, $random_gff_file) or die ("Unable to open file $random_gff_file: $!");
#	while(<$temp_data>){
#		print $out_data $_ unless($_ =~ /stddev=0.0;/);
#	}
#	close $temp_data;
#	close $out_data;
#	
#	#2 - merge + and - peaks
#	#sort by chr and by start coordinate first (required by bedtools)
#	system "sort -k1,1 -k4n,4 $random_gff_file > $random_sorted_gff_file";
#	if ( $? == -1 ){
#		print "sort: problem with output: $!\n";
#		exit -1;
#	}
#	#YOU NEED TO OPTIMISE THAT -d!!
#	system "$BEDTOOLS/mergeBed -n -d 35 -i $random_sorted_gff_file > $random_sorted_bed_file";
#	if ( $? == -1 ){
#		print "mergeBed: problem with output: $!\n";
#		exit -1;
#	}
#	
#	my $wc_out = `wc -l $random_sorted_bed_file`;
#    my @peaks = split(/ /, $wc_out);
#	#$peaks[0] -= 1 if ($peaks[0] > 0); #genetrack produces 1 line of header
#
#	return $peaks[0];
#}


sub run_gem{
	my ( $input_file, $gem_basename ) = @_;
	#system "java -Xmx10G -jar $GEM/gem.jar --g $CHROM_SIZE_PATH --d $GEM/Read_Distribution_ChIP-exo.txt --s $MAPPABILITY --expt $random_bam_file --f SAM --top 0 --outBED --out $random_gem_basename --genome $GENOME_PATH --k_min 2 --k_max 40";
	system "java -Xmx10G -jar $GEM/gem.jar --g $CHROM_SIZE_PATH --d $GEM/Read_Distribution_ChIP-exo.txt --s $MAPPABILITY --expt $input_file --f SAM --out $gem_basename";
	if ( $? == -1 ){
                print "GEM: problem with output: $!\n";
                exit -1;
    }
    #now capture number of rows in bedfile and fill an output file with them
    my $gem_peakfile = $gem_basename . "_GPS_events.txt";
    my $wc_out = `wc -l $gem_peakfile`;
    my @gem_peaks = split(/ /, $wc_out);
	return $gem_peaks[0];
}

#sub run_genetrack_single{
#
#	#python ${PCODE}/genetrack.py -v -s 5 -e 10 -F 6 -b ${FILE}
#	system "python $GENETRACK -s 5 -e 10 -F 6 $random_bed_file > $random_gff_file_s";
#        if ( $? == -1 ){
#                print "Genetrack: problem with output: $!\n";
#                exit -1;
#        }
#	my $wc_out = `wc -l $random_gff_file_s`;
#        my @genetrack_peaks = split(/ /, $wc_out);
#	$genetrack_peaks[0] -= 1 if ($genetrack_peaks[0] > 0); #genetrack produces 1 line of header
#
#	return $genetrack_peaks[0];
#}


#--nolambda            If True, MACS will use fixed background lambda as
#                        local lambda for every peak region. Normally, MACS
#                        calculates a dynamic local lambda to reflect the local
#                        bias due to potential chromatin structure.
#--llocal=LARGELOCAL   The large nearby region in basepairs to calculate
#                        dynamic lambda. This is used to capture the surround
#                        bias. DEFAULT: 10000.

sub run_macs14{
	my ( $input_file, $macs14_basename ) = @_;
	#/macs14 -t ${FILE} -n ${PDATA}/${ID} -s 40 --nomodel --shiftsize=13 --keep-dup=all -p 1e-8 -g 2540757438
	system "macs14 -t  $input_file -n $macs14_basename -g $MAPPABILITY -s 40 --keep-dup=all --nomodel --shiftsize 13 --verbose 1";	
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
	system "macs2 callpeak -t $input_file -n $macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 --verbose 1";	
	if ( $? == -1 ){
		print "MACS2: problem with output: $!\n";
                exit -1;
        }
	#now capture number of rows in bedfile and fill an outputfile with the
	my $macs2_bedfile = $macs2_basename . "_peaks.bed";
        my $wc_out = `wc -l $macs2_bedfile`;
        my @peaks = split(/ /, $wc_out);

	return $peaks[0];	
}





