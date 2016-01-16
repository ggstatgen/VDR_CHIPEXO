#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#needed to change the >1 dsdsddf entry to >1 paternal or >1 maternal
#from >6 dna:chromosome chromosome:GRCh37:6:1:171115067:1 to >6_paternal or >6_maternal


#notice you can use sed -e s/string1/string2/

my $infile;
GetOptions(
        'i=s'        =>\$infile
);
if(!$infile){
	print "USAGE: do_fasta_rename_headers.pl -i=<INFILE>\n";
    print "<INFILE> fasta file\n";
    exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	last if($_ =~ /^>MT/);
	
	my $line = $_;
	if($line =~ /^>(\d+)\s\w+/){
		$line = '>' . $1 . '_paternal';
	    #$line = '>' . $1 . '_maternal';
	}
	if($line =~ /^>(\w)\s\w+/){
		$line = '>' . $1 . '_paternal';
	    #$line = '>' . $1 . '_maternal';
	}
	print $line, "\n";
}
close $instream;