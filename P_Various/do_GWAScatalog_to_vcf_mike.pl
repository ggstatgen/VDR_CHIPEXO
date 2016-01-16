#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#INFO 14/11
#This converts the GWAS catalog .txt table prepared by Mike from the UCSC
#in a minimal vcf 
#The MINIMUM vcf understood by bedtools contains the following fields
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#So create these

#I need to intersect the genome wide (5E-8) snps coming from this with the LOOSE Ld blocks I generated here:
#/net/isi-mirror/1000_Genomes/POPULATION_PANELS/1kg_VCF/d_bed/plink_phase1_release_v3.20101123.EUR.bed.gz

#header
#bin    chrom   chromStart      chromEnd        name    pubMedID        author  pubDate journal title   trait   initSample      replSample      region  genes   riskAllele      riskAlFreq      pValue  pValueDesc      orOrBeta        ci95 platform        cnv

my $PROBABILITY = 0.00000005; 
my $MINSNPS = 5;

my $infile;
GetOptions(
        'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_GWAScatalog_to_vcf_mike.pl -i=<INFILE>\n";
     print "<INFILE> UCSC GWAS catalog TSV file\n";
     print "NOTE: I will only consider SNPs having p < $PROBABILITY and traits with at least $MINSNPS SNPs\n";
     print "CHANGE if needed\n";
     exit 1;
}

#print header
print "##fileformat=VCFv4.1\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");

#there will be duplicates (same all, different pvals)
#need a hash
my %textline;
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^#bin/); #header
	#vcf fields
	my $chrom; #1
	my $pos;   #3 chromend
	my $id;    #4
	my $missing = '-'; #REF     ALT     QUAL    
	my $filter; #10
	my $info;   #9
	my $pval; #17
	
	($info, $filter, $chrom, $pos, $id, $pval) = (split /\t/)[9,10,1,3,4,17];
	next if(!$chrom);
	next if($chrom eq '');
	
	$chrom = s/chr(.*)/$1/;
	
	#ONLY KEEP GENOME WIDE SIGNIFICANT SNPs
	next unless($pval);
	next if($pval eq 'NS'); #not sure what this means in the file
	next if($pval eq 'E');  # ''
	next if($pval > $PROBABILITY);
	$info = $info . '(p=' . $pval . ')';
	
	my $line =  $chrom . "\t" . $pos . "\t" . $id . "\t" . $missing . "\t" . $missing .  "\t" . $missing . "\t" . $filter . "\t" . $info . "\n";
	$textline{$line} = 1; 
}
close $instream;

#2nd pass
#build phenotype->snp-number map based on phenotype
my %phenotype_map;
foreach my $entry (keys %textline){
	my @fields = split (/\t/, $entry);
	if(!$phenotype_map{lc($fields[6])}){
		$phenotype_map{lc($fields[6])} = 1;
	}else{
		$phenotype_map{lc($fields[6])} += 1;		
	}
}

foreach my $item (keys %textline){
	my @fields = split (/\t/, $item);
	next if($phenotype_map{lc($fields[6])} < $MINSNPS); 
	print $item; 
}