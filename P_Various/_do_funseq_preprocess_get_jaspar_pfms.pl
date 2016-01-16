#OBSOLETE - used methodology in VDR POST-1REVIEW doc on google drive

#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#7/12/2015
#answering some of the reviewer questions
#Chris wants me to redo figure 3 using any other enriched Motifs

#I called again motifs using pscanchip, but this time instead of limiting myself to cpo3, I used as input a bed of all VDR-BVs and VDR-rBVs
#I select some PWMs which appear either locally, globally or centrally enriched according to pscanchip.
#I still don't have a good rationale on what to choose. I might simply pick 2 or 3 TFs that appeared strongly enriched after RXR:VDR in FIGURE s12

#to get the impact of the motif and alternative figure 3, I need to run funseq2. This means I need to add these motifs' pwms, and intervals, to funseq2

#This script converts the jaspar format for the pfms (positional frequency matrix) to the funseq format (positional frequency matrix normalised to one) and returns a file that needs to be appended to motif.PFM in the funseq2/data directory

#input format:
#>MA0074.1	RXRA::VDR
#A  [ 3  0  0  0  0  9  4  2  2  5  0  0  1  0  7 ]
#C  [ 0  0  0  0  9  0  2  4  0  0  0  0  0  9  1 ]
#G  [ 7 10  9  0  0  1  0  2  8  5 10  0  0  0  2 ]
#T  [ 0  0  1 10  1  0  4  2  0  0  0 10  9  1  0 ]

#output format (columns are A,C,G,T):
#>VDR_JASPAR MA0074.1 RXRA::VDR +
#R 0.300000 0.000000 0.700000 0.000000
#G 0.000000 0.000000 1.000000 0.000000
#G 0.000000 0.000000 0.900000 0.100000
#T 0.000000 0.000000 0.000000 1.000000
#C 0.000000 0.900000 0.000000 0.100000
#A 0.900000 0.000000 0.100000 0.000000
#N 0.400000 0.200000 0.000000 0.400000
#N 0.200000 0.400000 0.200000 0.200000
#R 0.200000 0.000000 0.800000 0.000000
#R 0.500000 0.000000 0.500000 0.000000
#G 0.000000 0.000000 1.000000 0.000000
#T 0.000000 0.000000 0.000000 1.000000
#T 0.100000 0.000000 0.000000 0.900000
#C 0.000000 0.900000 0.000000 0.100000
#A 0.700000 0.100000 0.200000 0.000000


#NOTE: I do not KNOW exactly how to attribute the letter in the first column. I will write a generic '.' character.

#db of all the input pfm matrices
my $INPUT_JASPAR = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/pfm_vertebrates.txt";


my $infile;
GetOptions(
        'i=s'        =>\$infile
);


$infile = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/infile_vdr_companion_pfms.txt";

if(!$infile){
	print "USAGE: do_pscanchip_out_intersect_vdrbvs.pl -i=<INFILE_IDs>\n";
    print "<INFILE_IDS> text file containing the jaspar IDs you want a normalised pfm for, one per line\n";
    exit 1;
}

my %jaspar_input_list;
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	$jaspar_input_list{$_} = 1;
}
close $instream;
my $found = scalar keys %jaspar_input_list;
print "Found $found unique Jaspar IDs in input file.\n";

###
#get all the jaspar matrices in an easy to access data structure
###

open ($instream,  q{<}, $INPUT_JASPAR) or die("Unable to open $INPUT_JASPAR : $!");
local $/;  # change the line separator to undef
my $jasparfile_contents = <$instream>;
close $instream;

my @matrices = split("\n>", $jasparfile_contents);
foreach my $item (@matrices){
	my %nt_A; my %nt_C; my %nt_G; my %nt_T;
	
	my ($identifier, $nt_A, $nt_C, $nt_G, $nt_T) = split("\n", $item);
	my ($id,$name) = split("\t", $identifier);
	$id =~ tr/>//d;
	
	#turn this into sub
	#clean up from '[]' characters
	$nt_A =~ tr/[]//d;
	$nt_C =~ tr/[]//d; 
	$nt_G =~ tr/[]//d;
	$nt_T =~ tr/[]//d;
	
	my @nt_A = split(" ", $nt_A);
	my @nt_C = split(" ", $nt_C);
	my @nt_G = split(" ", $nt_G);
	my @nt_T = split(" ", $nt_T);
	#remove the first
	shift @nt_A;
	shift @nt_C;
	shift @nt_G;
	shift @nt_T;
	
	my $CONSTANT1 = $nt_A[0] + $nt_C[0] + $nt_G[0] + $nt_T[0];
	my $CONSTANT2 = $nt_A[1] + $nt_C[1] + $nt_G[1] + $nt_T[1];
	
	print $item, "\n";
	print $CONSTANT1, "\t", $CONSTANT2,  "\n";
}


