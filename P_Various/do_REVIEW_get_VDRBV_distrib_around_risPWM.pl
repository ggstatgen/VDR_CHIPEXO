#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use List::Util qw( min max );

#07/2/2016
#I want to build a histogram showing the distribution of variant in the VDR:RXR motif and around it.
#I can use bedtools closest to assign every variant to its closest VDR:RXR motif and then spit the R file
#The problems with the interval in the ris are solved similarly to what I did in 'do_funseq_adapt_motiffile.pl'


#use the following
#bedtools closest -d -a Output_noDBRECUR.bed -b /net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_motif_scan_pscanchip/Pscanchip_occurrences_VDRRXR_jaspar_all_motifintervals.bed > closest.bed
#bedtools closest -D "ref" if you want it symmetrical

my $BIN_NUMBER = 1000;
my $BEDTOOLS = `which bedtools`; chomp $BEDTOOLS;
my $CHROMSIZES = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom_simple.sizes'; 

my $input_pscanchip_ris;
my $input_vdr_bv;
my $PWM_FILE;
my $identifier;
my $motif_name;
my $MIN_SCORE;
my %variant_binning;

GetOptions(
        'i_m=s'      =>\$input_pscanchip_ris,        
		'i_v'       =>\$input_vdr_bv,
        'pwm=s'		=>\$PWM_FILE,	
        'id=s'		=>\$identifier,
        'm=s'		=>\$motif_name,
        's=f'		=>\$MIN_SCORE        
);
$input_vdr_bv = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Output_noDBRECUR.bed";
$input_pscanchip_ris = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/RIS_LINKS/Pscanchip_hg19_bkgGM12865_Jaspar_VDRBVs_RXRA-VDR_MA0074.1_sites.ris";
$PWM_FILE = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Processed_PFMs_jaspar_FUNSEQ_INPUT.txt";
$motif_name = "RXRA-VDR";
$identifier = "MA0074.1";
$MIN_SCORE = '0.8';

my $USAGE = "\nUSAGE: $0 -i_m=<INFILE_PSCANCHIP> -i_v=<INFILE_VDRBV> -pwm=<ENCODE_PWM_FILE> -id=<ID> -m=<MOTIF_NAME> (opt)-s=<MINSCORE>\n" .
			"<INFILE_PSCANCHIP> ris file from PscanChip\n" .
			"<INFILE_VDRBV> bed file of VDR-BVs\n" . 
			"<ENCODE_PWM_FILE> text file with Jaspar PWMs in ENCODE format obtained with RSAT\n" .
			"<ID> string to use for the ID (e.g. Jaspar ID) of the PWM in the output file\n" . 
			"<MOTIF_NAME> string to use for the Motif PWM name in the output file\n"	.
			"optional <MINSCORE> lower threshold on score (eg 0.8) (default:none)\n"; 
			
unless($input_pscanchip_ris && $input_vdr_bv && $identifier && $motif_name && $PWM_FILE){
	print $USAGE;
	exit -1;
}
print STDERR "THRESHOLDING ON SCORE: $MIN_SCORE\n" if($MIN_SCORE);
my $this_motif_id = $motif_name . '_' . $identifier;

my($basename, $directory) = fileparse($input_pscanchip_ris);
$basename =~ s/(.*)\..*/$1/;
my $temp_pwm_bed        = $directory . $basename . '_' . $this_motif_id  . '_temp.bed';
my $temp_pwm_bed_sorted = $directory . $basename . '_' . $this_motif_id  . '_temp.sorted.bed';
my $data_closest        = $directory . $basename . '_' . $this_motif_id  . '_bedtools_closest.data';  
my $data_histogram      = $directory . $basename . '_' . $this_motif_id  . '_histogram.data';

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
system "$BEDTOOLS closest -D \"ref\" -a $input_vdr_bv -b $temp_pwm_bed_sorted -g $CHROMSIZES > $data_closest";
#unlink temp beds
unlink $temp_pwm_bed;
unlink $temp_pwm_bed_sorted;
#get histogram
system "cat $data_closest  | cut -f 7 | sort -k1,1n | uniq -c > $data_histogram";
unlink $data_closest;

#uniq -c has problems splitting the vars
exit;

bin_variants($data_histogram, \%variant_binning);



############### binning ############
sub bin_variants{
	my ($file, $hash) = @_;

	my $skipped = 0;
	my $valid = 0;
	my $total = 0;
	my %local_hash; #for the background the vcf is not unique, so fill hash. Do it also for the foreground, it won't hurt

# output from uniq -c
#      1 -6685283$
#      1 -6639883$
#      1 -6632412$
#      1 -6621090$
#      1 -6504254$
	#you need to find the maximum distance in absolute value
	my @distances_seen;
	open (my $instream,  q{<}, $file) or die("Unable to open $file : $!");
	while (<$instream>)	{
		chomp;
		my ($number, $distance) = (split / /)[7,8];
		$distance =~ s/-(\d+)/$1/;
		push(@distances_seen, $distance);
	}
	close $instream;
	my $max = max @distances_seen;
	
	open ($instream,  q{<}, $file) or die("Unable to open $file : $!");
	while (<$instream>)	{
		chomp;
		my ($number, $distance) = (split / /)[7,8];
		my $frac_distance = $distance / $max;
		my $bin_id = int( $frac_distance * $BIN_NUMBER )/ $BIN_NUMBER;
		$$hash{$bin_id} += $number;
	}	
	close $instream;	

	#my $bin_id = int($distance * $BIN_NUMBER) /$BIN_NUMBER;
	#$$hash{$bin_id} += $number;
	#my $bin_id = int( $variants_1kg{$item}{$AFGLOB} * $BIN_NUMBER ) / $BIN_NUMBER;

	return 1;
}
			
