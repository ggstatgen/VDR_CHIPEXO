#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#deepbind
#http://tools.genes.toronto.edu/deepbind/#

#This script expects a bed file as input (peaks or other). It will output a fasta file with the sequences you want use as Deepbind input. 

#if you run it on cpo3 or cpo10 etc, it will complain if any one of the fasta sequences is more than 1000 characters
#this script will trim sequences larger than 800bp by 100 left and 100 right

#Each sequence must be less than 1,000 bases long.
#Sequences of length 20-50 are recommended.


#syntax
#deepbind [--window-size n] [--average] [--echo] [--no-head] [--dump-info]
#         <id> [<id>...] [<sequence_file>]
#  Each <id> can be an actual ID (e.g. D00054.002), or it can be a file
#  name with extension .ids or .txt containing a single ID on each line.
#  
#  deepbind reads DNA/RNA sequences either from the <sequence_file>
#  or from stdin. If the first input character is '>', then FASTA input
#  format is assumed. Otherwise, the input is assumed to contain just raw
#  sequences, with one sequence on each line.
#  Each output line contains tab-separated scores, with
#  one score for each ID specified, in the order specified.
#  
#  --window-size n tells the model to use a length-n subsequence for each
#    per-position along the larger original input sequence; the default
#    is to use the length of the model's motif detectors (filters) x 1.5.
#  
#  --average causes the per-position scores to be averaged rather than maxed.
#  
#  --echo causes each output line echo the corresponding sequence;
#    ideal for pasting into a spreadsheet.
#  
#  --no-head omits the column names (IDs) from the first line of output.
#  
#  --dump-info causes information about each specified ID to be dumped to
#  stdout, e.g. the protein name, type, species, experiment, etc.


my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";
my $in_bed;
my $in_fasta;
my $MAX_INT = 500;

GetOptions(
        'i=s'      =>\$in_bed,
        'f=s'      =>\$in_fasta
);

#test 
#$in_bed = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/vdr_o3_peaks.bed";
#$in_fasta = "/net/isi-scratch/giuseppe/indexes/Hsap/hg19/hg19.fa";

if(!$in_bed){
	print "USAGE: do_deepbind_bed2fasta_sequences.pl -i=<IN_BED> -f=<IN_FASTA>\n";
    print "<IN_BED> bed file with input intervals\n";
    print "<IN_FASTA> reference FASTA to use for the output sequences\n";
    exit 1;
}
if(!$in_fasta){
	print "USAGE: do_deepbind_bed2fasta_sequences.pl -i=<IN_BED> -f=<IN_FASTA>\n";
    print "<IN_BED> bed file with input intervals\n";
    print "<IN_FASTA> reference FASTA to use for the output sequences\n";
    exit 1;
}
my($basename, $directory) = fileparse($in_bed);
$basename =~ s/(.*)\..*/$1/;

#prefilter the input sequences, and dump in a temp file
my $out_temp_bed =  $directory . $basename . '_temp.bed';
my $out_fasta = $directory . $basename . '.fasta';

print STDOUT "Trimming sequences larger than $MAX_INT to $MAX_INT..\n";

my %out_hash;
my $trim_count = 0;
open (my $instream,   q{<}, $in_bed) or die("Unable to open $in_bed : $!");
while(<$instream>){
	chomp;
	my @data = split("\t", $_);
	my $start = $data[1];
	my $stop = $data[2];
	
	if( ($stop - $start) >= $MAX_INT ){
		$trim_count++;
		my $length = $stop - $start;
		my $diff = int (($length - $MAX_INT)/2);
		$stop  = $stop - ($diff + 1);
		$start = $start + ($diff + 1);
	}
	
	my $out_line = $data[0] . "\t" . $start . "\t" . $stop;
	$out_hash{$out_line} = 1;
}
close $instream;
print STDOUT "Trimmed $trim_count input intervals.\n";

open (my $outstream,   q{>}, $out_temp_bed) or die("Unable to open $out_temp_bed : $!");
foreach my $item (keys %out_hash){
	print $outstream $item . "\n";
}
close $outstream;

#bedtools
#Usage:   bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf> -fo <fasta>  
system "$BEDTOOLS getfasta -fi $in_fasta -bed $out_temp_bed -fo $out_fasta";
print "Output save to $out_fasta\n";
print STDOUT "Finished.\n";
unlink $out_temp_bed;