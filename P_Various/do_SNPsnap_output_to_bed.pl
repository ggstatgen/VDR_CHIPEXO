#!/usr/bin/perl
use strict;
use warnings;
use IO::Zlib;

#3/7/2015
#Trying to get a frequency matched background for my foreground variables.
#I ran SNPsnap http://www.broadinstitute.org/mpg/snpsnap/
#with the 497 VDR-BV-REP input variants and got 10.000 random frequency matched sets of [497] SNPs

#How to couple this data with GAT?
#For now, get all the data and save in bed file
#Use background for GAT analyses

#a better way would be probably to intersect these snps with the background snps in peaks

#NOTE: tried intersecting the output of these (all SNPsnap snps) with SYMASYM background. The result is not frequency matched with the foreground.

my $INPUT_SNPSNAP = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/SNPsnap/SNPsnap_test_run/matched_snps.txt.gz";

#generate random index between 1 and 10,000
#my $minimum = 1;
#my $maximum = 10001;
#my $x = $minimum + int(rand($maximum));


tie *FILE,   'IO::Zlib', $INPUT_SNPSNAP, "rb";
while (<FILE>)	{
	chomp;
	next if($_ eq '');
	next if($_ =~ /^#/);
	next if($_ =~ /^Input_SNP/);
	
	my @line = split("\t", $_);
	
	shift @line; # first column is input data
	foreach my $item (@line){
		my ($chr,$pos) = split(":", $item);
		print $chr . "\t" . ($pos-1) . "\t" . $pos . "\n";	 		
	}
}
close FILE;

#my $data = $line[$x];
#my ($chr,$pos) = split(":", $data);
#print $chr . "\t" . ($pos-1) . "\t" . $pos . "\n";	 