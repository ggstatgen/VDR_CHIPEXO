#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;


#Script to prepare bed peak files for input to MEME
#Normally, you can get summits from MACS files, extend them +/- 150bp, get the underlying fasta files, pass them to MEME
#With my combined GEM/MACS approach, I no longer have summits, and have peaks of varying size.

#This file gets those peaks, gets the mid point, creates a bed output with only the coordinate of the midpoints, then calls bedtools to slop ALL of them by the same amount

my $infile_bed;
my $SLOP;
my $CHROM_SIZE_PATH = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes'; #change if you don't use hg19
GetOptions(
        'i=s'         =>\$infile_bed,
        's=i'         =>\$SLOP

);
if(!$infile_bed){
     print "USAGE: do_MEME_equalise_peak_width.pl -i=<BED_FILE> -s=<SLOP>\n";
     print "<BED_FILE> .bed file of peaks to process\n";
     print "<SLOP>  amount to extend inferred summit left and right\n";
     print "WARNING: slop refers to hg19, change if needed\n";
     exit 1;
}
if(!$SLOP){
     print "USAGE: do_MEME_equalise_peak_width.pl -i=<BED_FILE> -s=<SLOP>\n";
     print "<BED_FILE> .bed file of peaks to process\n";
     print "<SLOP>  amount to extend inferred summit left and right\n";
     print "WARNING: slop refers to hg19, change if needed\n";
     exit 1;
}
my $TOOL_PATH       = '/net/isi-scratch/giuseppe/tools';
my $BEDTOOLS        = $TOOL_PATH . '/bedtools2-2.20.1/bin';

my($basename, $directory) = fileparse($infile_bed);
$basename =~ s/(.*)\..*/$1/;
my $temp_out = $directory . $basename . '_summits.bed';

open (my $instream,     q{<}, $infile_bed) or die("Unable to open $infile_bed : $!");
open (my $tempstream,   q{>}, $temp_out) or die("Unable to open $temp_out : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /\#/);
	my ($chr, $start, $stop, $ID ) = (split /\t/)[0,1,2,3];
	next if (!$chr);
	
	my $summit = int( ($stop - $start)/2  );
	print $tempstream $chr, "\t", ($start + $summit - 1), "\t", ($start + $summit), "\t", $ID, "\n";
}
close $instream;
close $tempstream;

system "$BEDTOOLS/bedtools slop -g $CHROM_SIZE_PATH -b $SLOP -i $temp_out";
#| sort -k1,1V -k2,2g | $BEDTOOLS/bedtools merge -i stdin | $BEDTOOLS/bedtools coverage -abam $infile_bam -b stdin -counts > $full_idr_encodepeak_file_counts";
#unlink $temp_out;
