#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use Getopt::Long;

#GREAT does not like decimals in the signal field of the bed
#This script accepts as an input a normal bed from GEM or MACS
#outputs one compatible with GREAT


my $infile;
GetOptions(
	'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_bed_to_GREATcompatible_bed.pl -i=<FILE>\n";
     print "<FILE> input bed file\n";
     exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");

while(<$instream>){
	unless($_ =~ /^chr/){
		print $_;
		next;
	}
	my ($chrom, $chromStart, $chromEnd,$name,$score) = (split /\t/)[0, 1, 2,3,4];
	my $newline = $chrom . "\t" . $chromStart . "\t" . $chromEnd . "\t" . $name . "\t" . floor($score); 
	print $newline, "\n";
}
close $instream;
