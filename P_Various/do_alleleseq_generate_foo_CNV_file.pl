#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;


#24/11/2015
#Working on reviewers' output
#Need to run the Alleleseq pipeline on the VDR ChIPseq data 10855 and 10861
#I don't have full genome sequences and cannot run CNVnator to create RD values, unfortunately
#Therefore here I generate a CNV file from the SNV file with all rd set to 1.0
#output should be as follows
#chrm    snppos  rd
#1       52066   1
#1       695745  1
#1       742429  1

#input file structure
#1       23975   G       AG      CA      AA      HOMO
#1       38232   A       GG      GT      GG      HOMO
#1       38907   C       CT      TC      TT      HOMO
#1       41981   A       AG      GC      GG      HOMO
#1       46670   A       AA      AG      AG      HETERO

#just get 1 and 2 and add 1



my $infile_snp;
GetOptions(
	'i=s'  =>\$infile_snp,
);
if(!$infile_snp){
     print "USAGE: do_alleleseq_generate_foo_CNV_file.pl -i=<INFILE>\n";
     print "<INFILE> variants in the .snv file generated with do_alleleseq_vcf2snp_cluster.sh\n";
     exit 1;
}

my %snv_data;
open (my $instream,     q{<}, $infile_snp) or die("Unable to open $infile_snp : $!");
while(<$instream>){
	chomp;
	my ($chrm, $pos) = (split /\t/)[0,1];
	my $coord = $chrm . "\t" . $pos;
	$snv_data{$coord} = 1;
}
close $instream;

#print output
foreach my $item (sort keys %snv_data){
	print STDOUT $item . "\t" . '1' . "\n";
}

