#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#extract mapping reads of specified length
#input: desired length
#output: new bam file only containing those reads

my $infile;
my $outfile;
my $length;
GetOptions(
	"l=i" => \$length,
	'input=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_bam_extract_mappingreads_by_length.pl -i=<INFILE> -l=<LENGTH>\n";
     print "<INFILE> input bam file\n";
     print "<LENGTH> read length desired in the output\n";
     exit 1;
}
if(!$length){
     print "USAGE: gff_to_bed.pl -i=<INFILE>\n";
     print "<INFILE> input bam file\n";
     print "<LENGTH> read length desired in the output\n";
     exit 1;
}
#system "samtools view VDR_YRI_merged.sorted.bam | perl -ne '$rl_field = (split /\t/)[9]; print $_ if(length($rl_field) == 40);' | less -SN";