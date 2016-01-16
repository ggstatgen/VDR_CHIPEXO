#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Zlib;

#INFO 30/1/2015
#This converts the GWAS catalog .txt table found here
#http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz
#in a bed to be used with bedtools (I want to intersect it with my VDR-BVs to build forest plots for FIGURE 5 of the paper)


#you will need to save at least p-value, risk allele, OR, CI

#SAMPLE INPUT
#0 592	
#1 chr1
#2 1005805	
#3 1005806
#4 rs3934834
#5 19851299
#6 Johansson A
#7 2009-10-22
#8 Obesity (Silver Spring)
#9 Linkage and genome-wide association analysis of obesity-related phenotypes: association of weight with the MGAT1 gene.
#10 Body mass index
#11 1,079 South Tyrolian individuals, 790 Dutch founder individuals, 2,060 European ancestry individuals
#12 NA
#13 1p36.33
#14 NR
#15 rs3934834-G
#16 0.80
#17 6E-7
#18 (females + males)
#19 .11
#20 [NR] kg increase
#21 Illumina [318,237]
#22 N

#BED FIELDS
#chrom	start	end	name	score
#1	2	3	4-22(# separated)

my $infile;
my $b37_chromosomes;
GetOptions(
        'i=s'  =>\$infile,
        'b37'  =>\$b37_chromosomes,
);
if(!$infile){
     print "USAGE: do_GWAScatalog_UCSC_to_bed.pl -i=<INFILE> -b37\n";
     print "<INFILE> UCSC GWAS catalog zipped file\n";
     print "(optional)<b37> if you want b37 chromosome naming conventions in the output\n";
     exit 1;
}

tie *FILE,   'IO::Zlib', $infile, "rb";
while (<FILE>)	{ 
	chomp;
	next if($_ eq '');
	my @data = split("\t", $_);
	
	my $unknown = shift @data;
	my $chr = shift @data;
	$chr =~ s/chr(.*)/$1/ if($b37_chromosomes);
	
	my $start = shift @data;
	my $stop = shift @data;

	print STDOUT $chr . "\t" . $start . "\t" . $stop . "\t" . join('#', @data) . "\n";
}
close FILE;