#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#13/6/2013
#wrote this to clean up all the non "PASS" entries from a VCF file
#VCF 4.1
#VCF is a text file format (most likely stored in a compressed manner). 
#It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome.
#1- File meta-information is included after the ## string and must be key=value pairs.
#2 - The header line names the 8 fixed, mandatory columns after a #

my $infile; 
GetOptions(
	'i=s'  => \$infile,
);
if(!$infile){
     print "USAGE: clean_vcf.pl -i=<INFILE>\n";
     print "<INFILE> .vcf file to clean\n";
     exit 1;
}
my ($filename, $directory) = fileparse($infile);
$filename =~ s/(.*)\..*/$1/;
my $outfile =  $directory . "\/" . $filename . "_PASS.vcf";

open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
open (my $outstream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
while(<$instream>){
	#meta-info / header
	if($_ =~ /^\#/){
		print $outstream $_;
		next;
	}
	#data lines
	my $filter = (split /\t/)[6];
	if($filter =~ /PASS/){
		print $outstream $_;
	}else{
		next;
	}
}
close $instream;
close $outstream;
