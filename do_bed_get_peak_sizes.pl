#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use List::Util qw(min max);

#given an input bed, returns a tsv file with all peak lengths observed
#ID	LENGTH

#The idea is to use this after do_MACS2_merge_peaksets and before an R histogram analysis of the peak sizes
my $infile;
GetOptions(
        'i=s'         =>\$infile
);
if(!$infile){
     print "USAGE: do_bed_get_peaksizes.pl -i=<INFILE>\n";
     print "<INFILE> bed file with peak intervals\n";
     exit 1;
}

my($basename, $directory) = fileparse($infile);
my $ID;
$basename =~ s/(.*)\..*/$1/;

if($basename =~ /.*(GM\d{5}).*/){
	$ID = $1;
}else{
	$ID = $basename;
}

open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/); #header?
	next if($_ =~ /^track/); #header?
	my ($start, $stop) = (split /\t/)[1,2];
	#interesting to look by chromosome?
	my $length = $stop - $start;
	print $ID . "\t" . $length . "\n";
}
close $instream;
