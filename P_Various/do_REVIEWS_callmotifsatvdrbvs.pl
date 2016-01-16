#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#MODIFY to output fasta

#I want to try too see if there are VDR-BV falling in RXR:VDR motif we have not considered, because we were using the RXR:VDR motif in the CPo3 peaks and many VDR-BV are not in CPo3 peaks.

#This file accepts the Output.vcf from funseq2 and the size of the desired interval around the variant and outputs a bed file where each interval has a VDR-BV in the centre. This could be fed into PscanChip. Or alternatively into MEME-ChIP, though you would first need to convert to fasta


my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";
my $CHROM_SIZE_PATH = '/net/isi-scratch/giuseppe/indexes/chrominfo/g1k_v37_chrom.sizes';

my $input_file;
my $SLOP;
GetOptions(
        'i=s'        =>\$input_file,
        'slop=i'     =>\$SLOP
);

$SLOP = 50;
$input_file = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/Output_noDBRECUR.vcf";

my $USAGE = "USAGE: do_REVIEWS_callmotifsatvdrbvs.pl -i=<FUNSEQ_OUT> -slop=<INTERVAL>\n" .
			"<FUNSEQ_OUT> vcf output of funseq with vdr-bvs\n" .
			"<INTERVAL> size of the interval to consider around the variant\n";
			
if(!$input_file){
	print $USAGE;
    exit 1;
}
if(!$SLOP){
	print $USAGE;
    exit 1;
}

my($basename, $directory) = fileparse($input_file);
$basename =~ s/(.*)\..*/$1/;
my $temp_out  = $directory . 'TEMP_' . $basename . '.bed';
my $outfile  = $directory .  'REVIEWS_' . $basename . '_slop$INTERVAL.bed';

my %data;
open (my $instream,  q{<}, $input_file) or die("Unable to open $input_file : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/);
	next if($_ eq '');
	
	my ($chr, $pos) = (split /\t/)[0,1];
	my $coord = $chr . '-' . $pos;
	$data{$coord} = 1;
}
close $instream;

#create temp bed file and slop it
open (my $outstream,  q{>}, $temp_out) or die("Unable to open $temp_out : $!");
foreach my $item (keys %data){
	my ($chr, $coord) = split('-', $item);
	print $chr . "\t" . ($coord - 1) . "\t" . $coord . "\n";
}
close $outstream;

#slop it
#Usage:   bedtools slop [OPTIONS] -i <bed/gff/vcf> -g <genome> [-b <int> or (-l and -r)] 

#get fasta from paternal and maternal genomes
#system "$BEDTOOLS/bedtools slop -g $CHROM_SIZE_PATH -b $SLOP -i $out_bed | $BEDTOOLS/bedtools getfasta -fi $PATERNAL_GENOME_PATH -bed stdin -fo $out_fasta_pat";
system "$BEDTOOLS/bedtools slop -g $CHROM_SIZE_PATH -b $SLOP -i $temp_out > $outfile";


			