#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#quick script to turn the VDR peak file in a format compatible with the funseq annotation "ENCODE.annotation.gz"
#create this, then append to the existing encode annotation in funseq, sort, gzip

#format
#chr start stop TFP.VDR


my $infile;
GetOptions(
        'i=s'        =>\$infile
);
if(!$infile){
	print "USAGE: do_funseq_adapt_peakfile.pl -i=<INFILE>\n";
    print "<INFILE> bed file containing peaks\n";
    exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	my ($chr, $start, $stop) = (split /\t/)[0,1,2];
	
	print $chr, "\t", $start , "\t", $stop, "\t", 'TFP.VDR', "\n";
}
close $instream;