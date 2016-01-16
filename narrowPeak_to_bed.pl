#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#10/04
#Testing SPP on chip-exo files. Control is a dummy control with 100 reads selected at random from GM06986
#I want to see if the peaks called by SPP are comparable with others. But IGV won't accept SPP output peaks (narrowpeak)
#I convert them to bed
#23/04 also needed by all the other encode scripts

#from
#chr1    121485104       121485268       .       0       .       301.84984942945 -1      3.11182238802035        82
#to
#chr1    121485104       121485268       x

my $infile;
GetOptions(
	'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: narroPeak_to_bed.pl -i=<FILE>\n";
     print "<FILE> input narrowPeak file\n";
     exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
my $counter = 1;
while(<$instream>){
	my ($chrom, $chromStart, $chromEnd) = (split /\t/)[0, 1, 2];
	my $newline = $chrom . "\t" . $chromStart . "\t" . $chromEnd . "\t" . $counter; 
	print $newline, "\n";
	$counter++;
}
close $instream;
