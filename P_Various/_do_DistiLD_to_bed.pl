#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#INFO 15/11
#This converts the DistiLD .txt table found here
#http://distild.jensenlab.org/download.html
#into an annotated BED file for usage with GAT

#Andreas confirmed today that they don't use anymore the GWAS catalog (extending the snp position by 150KB left and right and turning it into a bed) but use this instead
#because it contains more precise LD blocks

#INPUT format
#0 PubMed ID of GWAS study
#1 Reference SNP (rs) number
#2 P-value
#3 Linkage disequilibrium (LD) block
#4 Ensembl genes in LD block
#5 Short description of GWAS study
#6 ICD10 codes (if applicable)

#eg
#0 21529783
#1 rs9556711
#2 8E-7
#3 chr13:97521524-98312811
#4 ENSG00000125249;ENSG00000125249;ENSG00000125249;ENSG00000165621;ENSG00000139793
#5 Alcoholism 12-month weekly alcohol consumption
#6 F10.2

#turn this into something of the form:
#CHROM	START	END	NAME	SCORE
#3-1	3-2	3-3	1(5)	2 

#OUTPUT 2 files:
#segment file for GAT --segments
#annotation name for GAT --descriptions


my $infile;
GetOptions(
        'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_DistiLD_to_bed.pl -i=<INFILE>\n";
     print "<INFILE> DistiLD TSV file\n";
     exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");

my ($filename, $directory) = fileparse($infile);
$filename =~ s/(.*)\..*/$1/;
my $outfile_segments = $directory . "\/" .  $filename . '_segments.bed';
my $outfile_descriptions = $directory . "\/" .  $filename . '_descriptions.txt';

#save into hash to remove duplicates
my %segments;
my %descriptions;
while(<$instream>){
	chomp;
	next if($_ eq '');
	#bed fields
	my $chr; my $start; my $end;
	my $ld_block; #3
	my $snp_id;    #1
	my $desc;  #5
	my $p_val; #2  

	($snp_id, $p_val, $ld_block, $desc) = (split /\t/)[1,2,3,5];
	if ($ld_block =~ /(chr\w+):(\d+)-(\d+)/){
		$chr = $1; $start = $2; $end = $3;
	}else{
		print STDERR "LD block: $ld_block in line $_ is in unrecognisable format.\n";
		#exit -1;
		next;
	}
	#build name field
	#my $name = $snp_id . ' (' . $desc . ')';
	#my $line = $chr . "\t" . $start . "\t" . $end . "\t" . $name . "\t" . $p_val . "\n";
	my $line = $chr . "\t" . $start . "\t" . $end . "\t" . $snp_id . "\t" . $p_val . "\n";
	$segments{$line} = 1;
	$descriptions{$snp_id} =  $desc;
}
close $instream;

open (my $segment_stream,  q{>}, $outfile_segments) or die("Unable to open $outfile_segments : $!");
foreach my $item (keys %segments){ print $segment_stream $item; }
close $segment_stream;

open (my $descriptions_stream,  q{>}, $outfile_descriptions) or die("Unable to open $outfile_descriptions : $!");
print $descriptions_stream "snp_id\tdescription\n";
foreach my $item (keys %descriptions){ print $descriptions_stream $item . "\t" . $descriptions{$item} . "\n"; }
close $descriptions_stream;