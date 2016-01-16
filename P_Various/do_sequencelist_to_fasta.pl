#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#from this
#abc
#cdu
#hgsfs
#sdfsf

#to this
#>seq1
#abc
#>seq2
#cdu
#>seq3
#hgsfs
#seq4
#sdfsf

#If the sequence is on the - strand, reverse complement it (WITH XXMOTIF IT'S ALREADY COMPLEMENTED)
#typical input
#funseq_VDRRXR_vdr_o10_peaks_meme_slop150.bed_Pvals.txt

my $SEQ_ROOT = ">seq_";
my $COUNTER = 1;

my $infile;
GetOptions(
        'i=s'        =>\$infile
);
if(!$infile){
	print "USAGE: do_sequencelist_to_fasta.pl -i=<INFILE>\n";
    print "<INFILE> output of XXmotif p-value tsv\n";
    exit 1;
}

open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	
	print $SEQ_ROOT . $COUNTER . "\n";
	print $_ . "\n";

	$COUNTER++;
}
close $instream;






sub reverse_complement_IUPAC {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}

