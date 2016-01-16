#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#standalone, or use in do_genetrack_cluster.sh

#27/02/2013
#PURPOSE 1
#I need to remove all lines of Genetrack output where singletons are found
#a singleton is indicated by genetrack with a std dev = 0.0
#eg 
#chr1    genetrack       .       568373  568383  6.24515535516   +       .       readcount=6;ID=P568378;stddev=2.35702260396;height=6.24515535516
#chr1    genetrack       .       568477  568487  6.0292452985    +       .       readcount=5;ID=P568482;stddev=0.0;height=6.0292452985
#input genetrack gff file
#output genetrack gff file without singletons

#PURPOSE 2
#I need to pair the peaks that Genetrack spits out
#peak(+)
#
#|-----|
#
#                        |-----|
#                         peak(-)
#
#             |
#             V
#             bedtools
#
#
#|------------------------------|
#    real peak 
#This probably means
#1 divide bed in two beds: 1)+ reads 2)- reads
#use bedtools windowbed?
#Using bedtools merge as recommended by the author in his email

#Tool:    bedtools merge (aka mergeBed)
#Version: v2.17.0
#Summary: Merges overlapping BED/GFF/VCF entries into a single interval.
#
#Usage:   bedtools merge [OPTIONS] -i <bed/gff/vcf>
#
#Options: 
#	-s	Force strandedness.  That is, only merge features
#		that are the same strand.
#		- By default, merging is done without respect to strand.
#
#	-n	Report the number of BED entries that were merged.
#		- Note: "1" is reported if no merging occurred.
#
#	-d	Maximum distance between features allowed for features
#		to be merged.
#		- Def. 0. That is, overlapping & book-ended features are merged.
#		- (INTEGER)
#
#	-nms	Report the names of the merged features separated by semicolons.
#
#	-scores	Report the scores of the merged features. Specify one of 
#		the following options for reporting scores:
#		  sum, min, max,
#		  mean, median, mode, antimode,
#		  collapse (i.e., print a semicolon-separated list),
#		- (INTEGER)
#
#Notes: 
#	(1) All output, regardless of input type (e.g., GFF or VCF)
#	    will in BED format with zero-based starts
#
#	(2) The input file (-i) file must be sorted by chrom, then start.

#sort the two beds as indicated here:
#http://bedtools.readthedocs.org/en/latest/content/tools/merge.html
#http://bedtools.readthedocs.org/en/latest/content/tools/sort.html

#ok splitting not necessary

my $infile;
my $tagdist; # tag distance for merging gff signals
my $BEDTOOLS = '/net/isi-scratch/giuseppe/tools/bedtools-2.17.0/bin';

GetOptions(
	'input=s'  =>\$infile,
	'd=s' =>\$tagdist,
);
if(!$infile){
     print "USAGE: genetrack_postprocessing.pl -input=<INFILE> -d=<TAGDIST>\n";
     print "<INFILE> input fastq file\n";
     print "<TAGDIST> tag distance for mergeBed (e.g. 35)\n";
     exit 1;
}
if(!$tagdist){
     print "USAGE: genetrack_postprocessing.pl -input=<INFILE> -d=<TAGDIST>\n";
     print "<INFILE> input fastq file\n";
     print "<TAGDIST> tag distance for mergeBed (e.g. 35)\n";
     exit 1;
}
my $file_nosingletons = $infile;
$file_nosingletons =~ s/(.*)\..*/$1/;
my $file_nosingletons_sorted = $file_nosingletons . "_sorted_nosgtons.gff";
my $file_nosingletons_merged = $file_nosingletons . ".bed";
$file_nosingletons = $file_nosingletons . "_nosgtons.gff";

open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
open (my $outstream,  q{>}, $file_nosingletons) or die("Unable to open $file_nosingletons : $!");
while(<$instream>){
	print $outstream $_ unless($_ =~ /stddev=0.0;/);
	my $strand = (split /\t/)[6];
}
close $instream;
close $outstream;

#system "$BEDTOOLS/sortBed -i $file_nosingletons > $file_nosingletons_sorted";
system "sort -k1,1 -k4n,4 $file_nosingletons > $file_nosingletons_sorted";
if ( $? == -1 ){
	print "sort: problem with output: $!\n";
	exit -1;
}
#YOU NEED TO OPTIMISE THAT $tagdist!!
system "$BEDTOOLS/mergeBed -n -d $tagdist -i $file_nosingletons_sorted > $file_nosingletons_merged";
if ( $? == -1 ){
	print "mergeBed: problem with output: $!\n";
	exit -1;
}
unlink $file_nosingletons;
unlink $file_nosingletons_sorted;

