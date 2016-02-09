#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#07/2/2016
#I want to build a histogram showing the distribution of variant in the VDR:RXR motif and around it.
#I can use bedtools closest to assign every variant to its closest VDR:RXR motif and then spit the R file
#The problems with the interval in the ris are solved similarly to what I did in 'do_funseq_adapt_motiffile.pl'

my $BEDTOOLS = `which bedtools`; chomp $BEDTOOLS;
my $RSCRIPT = `which RRscript`; chomp $RSCRIPT;

my $CHROMSIZES = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom_simple.sizes'; 
my $IN_VDRBV_BED = '/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Output_noDBRECUR.bed';
my $IN_VDRrBV_BED = '/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Output_VDR_rBVs_hg19.bed';
my $PWM_FILE = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Processed_PFMs_jaspar_FUNSEQ_INPUT.txt";

my $input_pscanchip_ris;
my $input_variants;
my $INPUT_VAR_BINARY;
my $identifier;
my $motif_name;
my $THRS_DIST;
my $PLOT_EXT;
my $MIN_SCORE;

GetOptions(
        'm=s'		=>\$input_pscanchip_ris,        
		'v=s'		=>\$INPUT_VAR_BINARY,
        't=i'		=>\$THRS_DIST,   
        'p=s'		=>\$PLOT_EXT,
        's=f'		=>\$MIN_SCORE            
);

my $USAGE = "\nUSAGE: $0 -m=<INFILE_PSCANCHIP> -v=<BV|rBV> -t=<THRS> -p=<pdf|svg|jpg|png> (opt)-s=<MINSCORE>\n" .
			"<INFILE_PSCANCHIP> ris file from PscanChip\n" .
			"<BV|rBV> if BV, all VDRBV will be used; if rBV, only VDR-rBV will be used\n" .
			"<THRS> Distance threshold cutoff for plot\n" .
			"<PLOT> type of R output plot desired\n" .
			"optional <MINSCORE> lower threshold on score (eg 0.8) (default:none)\n"; 
			
unless($input_pscanchip_ris && $INPUT_VAR_BINARY && $identifier && $motif_name && $THRS_DIST && $PLOT_EXT){
	print $USAGE;
	exit -1;
}
#binary arguments
if($INPUT_VAR_BINARY eq 'BV'){
	$input_variants = $IN_VDRBV_BED;
}elsif($INPUT_VAR_BINARY eq 'rBV'){
	$input_variants = $IN_VDRrBV_BED;
}else{
	print STDERR "ERROR: field -v not recognised: $INPUT_VAR_BINARY. Aborting.\n";
	exit -1;
}
unless($PLOT_EXT eq 'jpg' || $PLOT_EXT eq 'png' || $PLOT_EXT eq 'svg' || $PLOT_EXT eq 'pdf'){
	print STDERR "ERROR: field -p not recognised: $PLOT_EXT. Aborting.\n";
	exit -1;
}
print STDERR "THRESHOLDING ON SCORE: $MIN_SCORE\n" if($MIN_SCORE);

$INPUT_VAR_BINARY = 'VDR-' . $INPUT_VAR_BINARY;
my($basename, $directory) = fileparse($input_pscanchip_ris);
$basename =~ s/(.*)\..*/$1/;
#get motif name and id from the filename
#eg Pscanchip_hg19_bkgGM12865_Jaspar_VDRBVs_RXRA-VDR_MA0074.1_sites
if($basename =~ /BVs_(.*)_(.*)_sites/){
	$motif_name = $1;
	$identifier = $2;
}else{
	print STDERR "ERROR: Unable to recognise the input motif file name: $basename. Aborting.\n";
	exit -1;
}

my $this_motif_id = $motif_name . '_' . $identifier;
my $temp_pwm_bed            = $directory . $basename . '_' . $this_motif_id  . '_temp.bed';
my $temp_pwm_bed_sorted     = $directory . $basename . '_' . $this_motif_id  . '_temp.sorted.bed';
my $data_closest            = $directory . $basename . '_' . $this_motif_id  . '_bedtools_closest.data';  
my $data_counts             = $directory . 'R_' . $this_motif_id  . '_counts.Rdata';
my $data_histogram          = $directory . 'R_' . $this_motif_id  . '_histogram.Rdata';

my $Rscript_counts           = $directory . 'R_' .                     $this_motif_id  . '_counts.R';
my $Rscript_hist             = $directory . 'R_' . $this_motif_id  . '_hist.R';

my $Rscript_counts_plot_all  = $directory . 'R_' . $INPUT_VAR_BINARY . '_' . $this_motif_id . '_counts_all.'                       . $PLOT_EXT;
my $Rscript_counts_plot_sub  = $directory . 'R_' . $INPUT_VAR_BINARY . '_' . $this_motif_id . '_counts_dthrs_' . $THRS_DIST .  '.' . $PLOT_EXT;
my $Rscript_hist_plot_all    = $directory . 'R_' . $INPUT_VAR_BINARY . '_' . $this_motif_id . '_hist_all.'                         . $PLOT_EXT;
my $Rscript_hist_plot_sub    = $directory . 'R_' . $INPUT_VAR_BINARY . '_' . $this_motif_id . '_hist_dthrs_'   . $THRS_DIST .  '.' . $PLOT_EXT;

#######
#1 get the REAL motif length from the encode representations of the motif
#######
my %motif; my $A = 1; my $C = 2; my $G = 3; my $T = 4;
my $prev_name; my @info; my $temp;
#slurp pwm file
open (my $instream,  q{<}, $PWM_FILE) or die("Unable to open $PWM_FILE : $!");
        while(<$instream>){
                chomp $_;
                if(/^>/){
                        $prev_name = (split/>|\s+/,$_)[1];
                }else{
                        @info = split/\s+/,$_;
                        if(not exists $motif{$prev_name}){
                                $motif{$prev_name}->[0] = {(A=>$info[$A], T=>$info[$T], C=>$info[$C], G=>$info[$G])};
                        }else{
                                $temp = $motif{$prev_name};
                                $motif{$prev_name}->[scalar(@$temp)] = {(A=>$info[$A], T=>$info[$T], C=>$info[$C], G=>$info[$G])};
                        }
                }
        }
close $instream;

#heterodimers are saved by Jaspar as monomer::monomer
#I replaced the :: with a '-' in the input file name because it's not recognised by the SGE submission
#change again here:
$motif_name =~ s/\-/\:\:/;
#get the length of the motif analyzed in this iteration
my $full_motif_id = $motif_name . '_' . $identifier;
my $ref = $motif{$full_motif_id};
my $MOTIF_LENGTH = scalar(@$ref); 
print STDERR "The length of the motif: $full_motif_id according to the JASPAR Pwm is $MOTIF_LENGTH\n";

my %pwm_interval_bed;
open ($instream,  q{<}, $input_pscanchip_ris) or die("Unable to open $input_pscanchip_ris : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^CHR/);
	
	my ($chr,$motif_start,$motif_end,$motif_strand,$score,$site) = (split /\t/)[0,4,5,8,9,10];
	#next if (!$chr);
	if(  $MIN_SCORE && ($score < $MIN_SCORE) ){
		next;
	}
	my $pscanchip_interval_length = ($motif_end - $motif_start);

	if($pscanchip_interval_length <  $MOTIF_LENGTH){
		$motif_end += 1;
	}elsif($pscanchip_interval_length == $MOTIF_LENGTH){	
		;
	}else{
		print STDERR "ERROR: the pscanchip motif length is LARGER than the Jaspar length. Verify.\n";
		exit -1;	
	}
	my $bed_line = $chr . "\t" . $motif_start . "\t" . $motif_end;
	$pwm_interval_bed{$bed_line} = 1;
}
close $instream;

open (my $outstream,  q{>}, $temp_pwm_bed) or die("Unable to open $temp_pwm_bed : $!");
foreach my $item (keys %pwm_interval_bed){ print $outstream $item, "\n"; }
close $outstream;

#sort bed
system "cat $temp_pwm_bed | sort -k1,1V -k2,2n | uniq > $temp_pwm_bed_sorted";
#get closeness data
system "$BEDTOOLS closest -D \"ref\" -a $input_variants -b $temp_pwm_bed_sorted -g $CHROMSIZES > $data_closest";
#from here I want two files
#one, pure counts, to build histograms in R
#two, x,y pairs, to build point plots in R

#get pure counts
system "cat $data_closest  | cut -f 7 > $data_histogram";
#get x,y pairs, where x=distance and y=count at that distance
system "cat $data_closest  | cut -f 7 | sort -k1,1n | uniq -c | sed \'s/^ *//g;s/ /\t/\'  > $data_counts";

unlink $temp_pwm_bed;
unlink $temp_pwm_bed_sorted;
unlink $data_closest;

#line plot
open ($outstream,  q{>}, $Rscript_counts) or die("Unable to open $Rscript_counts : $!");
print $outstream "data <- read.table(\"$data_counts\",sep=\"\\t\")" . "\n";
print $outstream "$PLOT_EXT(file=\"$Rscript_counts_plot_all\")" . "\n";
print $outstream "plot(data\$V2,data\$V1, type=\"p\", cex=.5, xlab=\"Distance from meta-motif PWM (bp)\", ylab=\"\#Observations\", main=\"$INPUT_VAR_BINARY profile around $full_motif_id\")" . "\n";
print $outstream "dev.off()" . "\n";
print $outstream "sub_data <- subset(data, V2 >= -$THRS_DIST & V2 <= $THRS_DIST, select=c(V1,V2))" . "\n";
print $outstream "$PLOT_EXT(file=\"$Rscript_counts_plot_sub\")" . "\n";
print $outstream "plot(sub_data\$V2,sub_data\$V1, type=\"p\", cex=.3, xlab=\"Distance from meta-motif PWM (bp)\", ylab=\"\#Observations\", main=\"$INPUT_VAR_BINARY profile around $full_motif_id\")" . "\n";
print $outstream "dev.off()" . "\n";
close $outstream;

system "$RSCRIPT $Rscript_counts";

open ($outstream,  q{>}, $Rscript_hist) or die("Unable to open $Rscript_hist : $!");
print $outstream "data <- read.table(\"$data_histogram\",sep=\"\\t\")" . "\n";
print $outstream "$PLOT_EXT(file=\"$Rscript_hist_plot_all\")" . "\n";
print $outstream "hist(data\$V1, 100, xlab=\"Distance from meta-motif PWM (bp)\", ylab=\"\Frequency\", main=\"$INPUT_VAR_BINARY distribution around $full_motif_id\")" . "\n";
print $outstream "dev.off()" . "\n";
print $outstream "sub_data <- subset(data, V1 >= -$THRS_DIST & V1 <= $THRS_DIST)" . "\n";
print $outstream "$PLOT_EXT(file=\"$Rscript_hist_plot_sub\")" . "\n";
print $outstream "hist(sub_data\$V1, 100, xlab=\"Distance from meta-motif PWM (bp)\", ylab=\"\Frequency\", main=\"$INPUT_VAR_BINARY distribution around $full_motif_id\")" . "\n";
print $outstream "dev.off()" . "\n";
close $outstream;

system "$RSCRIPT $Rscript_hist";

unlink $data_counts;
unlink $data_histogram;
unlink $Rscript_counts;
unlink $Rscript_hist;




#bin_variants($data_histogram, \%variant_binning);
#print "BIN\tCOUNT\n";
#foreach my $bin (sort {$a<=>$b} keys %variant_binning){
#	print $bin, "\t", $variant_binning{$bin}, "\n";
#}
################ binning ############
#sub bin_variants{
#	my $BIN_NUMBER; #temp
#	my ($file, $hash) = @_;
#	my %local_hash; #for the background the vcf is not unique, so fill hash. Do it also for the foreground, it won't hurt
#
## output from uniq -c
##      1 -6685283$
##      1 -6639883$
##      1 -6632412$
##      1 -6621090$
##      1 -6504254$
#	#you need to find the maximum distance in absolute value
#	my @distances_seen;
#	open (my $instream,  q{<}, $file) or die("Unable to open $file : $!");
#	while (<$instream>)	{
#		chomp;
#		my ($number, $distance) = (split /\t/)[0,1];
#		$distance =~ s/-(\d+)/$1/;
#		push(@distances_seen, $distance);
#	}
#	close $instream;
#	my $max = max @distances_seen;
#	#print "Max distance in absolute value: $max\n";
#	
#	open ($instream,  q{<}, $file) or die("Unable to open $file : $!");
#	while (<$instream>)	{
#		chomp;
#		my ($number, $distance) = (split /\t/)[0,1];
#		my $frac_distance = $distance / $max;
#		my $bin_id = int( $frac_distance * $BIN_NUMBER )/ $BIN_NUMBER;
#		$$hash{$bin_id} += $number;
#	}	
#	close $instream;	
#	return 1;
#}
