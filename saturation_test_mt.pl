#!/usr/bin/perl
use strict;
use warnings;
use threads;
use Getopt::Long;
use File::Basename;
use List::Util qw(sum);

#14/02/2013
#I want to see if there is a minimum number of reads that guarantee maximum numbers of peaks found
#In other words I want to plot the number of peaks found against a set of reads of increasing size and see
#do we observe saturation?
#This is important because we have many files with low read count. 
#What's the minimum number of reads to get similar number of peaks?
#Do we get saturation?

#This script should
#1 obtain a subset of randomly selected reads from a fastq file
#2 map them to human genome
#3 find peaks and obtain the total number of peaks
#cycle over again, with a higher subset of reads
#stop when all the file has been used.

#INPUTS
#1 fastq file
#2 step or frequency (0.01? 0.1)

#OUTPUTS
#TSV FILE  <number of reads> <uniquely mapped reads> <number of peaks> <number of peaks/total peaks for full file>

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
#temp data
#---------
my $infile_sam = $datadir . "\/" . $basename . '.sam';
my $infile_bam = $datadir . "\/" . $basename . '.bam';
my $infile_bai = $datadir . "\/" . $basename . '.bai';
my $infile_bed = $datadir . "\/" . $basename . '.bed';

my $outfile = $datadir . "\/" . $basename . '.result';

my $random_fastq_file       = $datadir . "\/" . $input_id  . '_random.fq'; #the python fastq sumbsampler writes here / bowtie reads
my $random_sam_file         = $datadir . "\/" . $input_id  . '_random.sam'; #bowtie writes here / picard reads
my $random_bam_file         = $datadir . "\/" . $input_id  . '_random.bam'; #picard writes here /macs2 and GEM read
my $random_bai_file         = $datadir . "\/" . $input_id  . '_random.bai'; #picard writes here
my $random_bed_file         = $datadir . "\/" . $input_id  . '_random.bed'; #bamtobed writes here / genetrack reads
my $random_macs2_basename   = $datadir . "\/" . $input_id  . '_random_macs2'; #macs2 outputs in these
my $random_macs14_basename  = $datadir . "\/" . $input_id  . '_random_macs14'; #macs14 outputs in these
my $random_gem_basename     = $datadir . "\/" . $input_id  . '_random_gem'; #gem outputs in these
#my $random_gff_file_s         = $input_id  . '_random.gff'; #genetrack outputs here

#unlink files in case something is left from previous run
unlink $random_fastq_file; unlink $random_sam_file; unlink $random_bam_file; unlink $random_bai_file; unlink $random_bed_file;
unlink $infile_sam; unlink $infile_bam; unlink $infile_bed;
#unlink $random_gff_file_s;

#my $wc_out = `wc -l $infile`;
#my @wc_out = split(/ /, $wc_out);
#my $size = $wc_out[0]/4;#number of lines in single fastq entry
#print '-----Initial .fastq has ' . $size . ' entries' . "\n"; 

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
#bam to bed (only needed by genetrack threads - improve )
#------------
#convert bam to bed
#${PCODE}/bamToBed -i ${FILE} > ${PDATA}/${BASEFILE}.bed;
system "$BEDTOOLS/bamToBed -i $infile_bam > $infile_bed";
if ( $? == -1 ){
	print "bamToBed: problem with output: $!\n";
    exit -1;
}
#---------------
#peak calling
#---------------
#genetrack threads, one per chromosome
for my $i (0 .. $#HUMAN_CHRS) {
	push(@full_threads, threads->create(\&run_genetrack, $HUMAN_CHRS[$i], $infile_bed));
}
foreach my $thread (@full_threads){
	my $peak = $thread->join;
	push(@full_peaks, $peak);	
}
my $total_genetrack_peaks = sum(@full_peaks);
print "Genetrack: number of total peaks in initial file: " . $total_genetrack_peaks . "\n";

#remove all temp files
unlink $infile_sam;
unlink $infile_bam;
unlink $infile_bai;
unlink $infile_bed;
my $gt_outputs_gff = $datadir . "\/" . "$input_id\*chr\*.gff";
my $gt_outputs_bed = $datadir . "\/" . "$input_id\*chr\*.bed";
system "rm $gt_outputs_gff"; #rm genetrack outputs
system "rm $gt_outputs_bed"; #rm genetrack outputs
#########################
#randomisation main loop
#########################
#outstream header
open (my $outstream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
#print $outstream "SAMPLE\tREADS\tPEAKS_MACS\tPEAKS_GENETRACK\tPEAKS_GEM\n";
#print $outstream "SAMPLE\tREADS\tPEAKS_MACS14\tPEAKS_MACS2\tPEAKS_GENETRACK\tPEAKS_GEM\n";
#print $outstream "SAMPLE\tREADS\tPEAKS_GENETRACK\n";
print $outstream "SAMPLE\tREADS_FQ\tMAPPED_READS_BAM\tPEAKS_GENETRACK\tPEAK_RATIO_GENETRACK\n";
my $counter = 1;
while ($INIT <= $MAX){
	my $wc_out; my @wc_out;
	my @gt_threads; my @gt_peaks;
	#----------------
	#random subset selection
	#----------------
	print "Creating random .fastq file, (fraction: $INIT)..\n";
	#usage: subsample.py <fraction> <input file> <output file>
	system "$SCRIPT_PATH/do_subsample_single.py $INIT $infile $random_fastq_file";

	$wc_out = `wc -l $random_fastq_file`;
	@wc_out = split(/ /, $wc_out);
	my $number_of_fq_reads = $wc_out[0]/4;	
    print '------Random .fastq sample number ' . $counter . ' has ' . $number_of_fq_reads . ' entries' . "\n";
	#--------------
	#align
	#-------------
	print "Aligning random .fastq file..\n";
	system "$BOWTIE -q -S --best --strata -m 1 -p 12 --chunkmbs 1024 $INDEX_PATH $random_fastq_file $random_sam_file";
	if ( $? == -1 ){
		print "BOWTIE: problem with output: $!\n";
		exit -1;
	}
	#-------------
	#sam to bam
	#------------
    #java -Xmx4g -Djava.io.tmpdir=${TMP} -jar ${PCODE_PICARD}/SortSam.jar SORT_ORDER=coordinate INPUT=${FILE} OUTPUT=${PDATA}/${FILE_ID}.bam CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"
	system "java -Xmx4g -Djava.io.tmpdir=$TMP_DIR -jar $PICARD/SortSam.jar SORT_ORDER=coordinate INPUT=$random_sam_file OUTPUT=$random_bam_file CREATE_INDEX=true VERBOSITY=WARNING VALIDATION_STRINGENCY=LENIENT";
	my $number_of_mapped_reads = `$SAMTOOLS/samtools view -F 4 $random_bam_file | wc -l`; #uniquely mapped reads
	chomp $number_of_mapped_reads;
	#------------
	#bam to bed (only needed by genetrack threads - improve )
	#------------
	#${PCODE}/bamToBed -i ${FILE} > ${PDATA}/${BASEFILE}.bed;
	system "$BEDTOOLS/bamToBed -i $random_bam_file > $random_bed_file";
	if ( $? == -1 ){
                print "bamToBed: problem with output: $!\n";
                exit -1;
	}
	#---------------
	#peak calling
	#---------------
	#my $tmacs14 = threads->create( \&run_macs14 );
	#my $tmacs2  = threads->create( \&run_macs2 );
	#my $tgps    = threads->create( \&run_gem );
	
	#my $t_gt  = threads->create( \&run_genetrack_single );
	#my $gtsingle_peaks  = $t_gt->join;
	#genetrack threads, one per chromosome
	for my $i (0 .. $#HUMAN_CHRS) {
		push(@gt_threads, threads->create(\&run_genetrack, $HUMAN_CHRS[$i], $random_bed_file));
	}
	#push (@gt_peaks, $_->join) foreach @gt_threads;
	foreach my $thread (@gt_threads){
		my $peak = $thread->join;
		push(@gt_peaks, $peak);	
	}
	#my $macs14_peaks = $tmacs14->join;
	#my $macs2_peaks  = $tmacs2->join;
	my $genetrack_peaks = sum(@gt_peaks);
	#my $gem_peaks = $tgps->join;

	print $outstream  $counter . "\t" . $number_of_fq_reads . "\t" .  $number_of_mapped_reads . "\t" .  $genetrack_peaks .  "\t" . ($genetrack_peaks/$total_genetrack_peaks) .  "\n";
	#print $counter . "\t" . $number_of_reads . "\t" . $macs14_peaks . "\t" . $macs2_peaks . "\t" . $genetrack_peaks . "\n";
	#print $outstream  $counter . "\t" . $number_of_reads . "\t" . $macs14_peaks . "\t" . $macs2_peaks . "\t" . $genetrack_peaks .  "\t" . $gem_peaks . "\n";

	$INIT += $STEP;
	$counter += 1;
	unlink $random_fastq_file;
	unlink $random_sam_file;
	unlink $random_bam_file;
	unlink $random_bai_file;
	unlink $random_bed_file;

   	#unlink $random_gff_file_s;
	system "rm $gt_outputs_gff"; #rm genetrack outputs
	system "rm $gt_outputs_bed"; #rm genetrack outputs
	system "rm $random_macs2_basename*.*"; #rm macs2 output
	system "rm $random_macs14_basename*.*"; #rm macs14 output
	system "rm -rf $random_gem_basename*"; #remove gem output
}

print $outstream  $counter . "\t" . $total_fq_reads . "\t" .  $total_mapped_reads . "\t" .  $total_genetrack_peaks .  "\t" . ($total_genetrack_peaks/$total_genetrack_peaks) .  "\n";
close $outstream;


#this will need
#its own chromosome
#its own random file
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
	
	my $wc_out = `wc -l $random_sorted_bed_file`;
    my @peaks = split(/ /, $wc_out);
	#$peaks[0] -= 1 if ($peaks[0] > 0); #genetrack produces 1 line of header

	return $peaks[0];
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
	#/macs14 -t ${FILE} -n ${PDATA}/${ID} -s 40 --nomodel --shiftsize=13 --keep-dup=all -p 1e-8 -g 2540757438
	system "macs14 -t $random_bam_file -n $random_macs14_basename -g $MAPPABILITY -s 40 --keep-dup=all --nomodel --shiftsize 13 --verbose 1 --nolambda --llocal 0";	
	if ( $? == -1 ){
		print "MACS14: problem with output: $!\n";
                exit -1;
        }
	#now capture number of rows in bedfile and fill an outputfile with the
	my $macs14_bedfile = $random_macs14_basename . "_peaks.bed";
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
	#${PCODE}/macs2 callpeak -t ${FILE} -n ${PDATA}/${ID} -g 2540757438 --keep-dup=all -q 0.001 --nomodel --extsize 26 --call-summits
	system "macs2 callpeak -t $random_bam_file -n $random_macs2_basename -g $MAPPABILITY --keep-dup=all --nomodel --extsize 26 --verbose 1 --nolambda --llocal 0";	
	if ( $? == -1 ){
		print "MACS2: problem with output: $!\n";
                exit -1;
        }
	#now capture number of rows in bedfile and fill an outputfile with the
	my $macs2_bedfile = $random_macs2_basename . "_peaks.bed";
        my $wc_out = `wc -l $macs2_bedfile`;
        my @peaks = split(/ /, $wc_out);

	return $peaks[0];	
}




sub run_gem{
	#system "java -Xmx10G -jar $GEM/gem.jar --g $CHROM_SIZE_PATH --d $GEM/Read_Distribution_ChIP-exo.txt --s $MAPPABILITY --expt $random_bam_file --f SAM --outBED --out $random_gem_basename --genome $GENOME_PATH --k_min 2 --k_max 40";
	system "java -Xmx10G -jar $GEM/gem.jar --g $CHROM_SIZE_PATH --d $GEM/Read_Distribution_ChIP-exo.txt --s $MAPPABILITY --expt $random_bam_file --f SAM --t 10 --out $random_gem_basename";
	if ( $? == -1 ){
                print "GEM: problem with output: $!\n";
                exit -1;
        }
        #now capture number of rows in bedfile and fill an outputfile with them
        my $gem_peakfile = $random_gem_basename . "_GPS_events.txt";
        my $wc_out = `wc -l $gem_peakfile`;
        my @gem_peaks = split(/ /, $wc_out);
	return $gem_peaks[0];
}
