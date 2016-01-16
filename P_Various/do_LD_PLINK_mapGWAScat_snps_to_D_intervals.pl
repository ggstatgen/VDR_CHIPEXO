#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use IO::Zlib;

#Figure 6
#map gwascat snps to D' ld blocks

#Given the 1000g LD blocks for CEU and YRI obtained with Plink  and stored here
#/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_CEU_LD_b37.bed.gz

#I want to create a bed where I replace, for each GWAS snp reaching genome-wide significance, its LD block.
#crucially, if I find a GWAS snp does not intersect any LD block, the output bed will feature a 1bp interval with that GWAS snp's coordinates

#LD INTERVALS
my $INPUT_LD_CEU = "/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_CEU_LD_b37.bed.gz";
#my $INPUT_LD_YRI = "/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_YRI_LD_b37.bed.gz";

#DISEASE INFO
#for both, only snps with p < 5e-8; only traits with 5+ snps at least
#18/2/2015 now using up to date UCSC GWAS catgwasCatalog.b37.maxp_5E-8.vcf.gz
my $INPUT_GWAS = "/net/isi-scratch/giuseppe/indexes/GWAS_Gwascat/UCSC_gwas_catalog/gwasCatalog.hg19.Beecham92-MS.maxp_5E-8.vcf.gz";
my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";

print STDERR "do_LD_gwassnp2block.pl:\n";
print STDERR "Using the following files:\n";
print STDERR "gwas catalog vcf: $INPUT_GWAS\n";
print STDERR "LD blocks: $INPUT_LD_CEU\n";
print STDERR "EDIT if required\n";

#output directory will be where the disease data is 
my ($filename_gwas,  $out_dir_gwas) = fileparse($INPUT_GWAS);
print STDERR "Saving all output under $out_dir_gwas..\n";

#temp: all gwas snps not falling in LD blocs
my $temp_vcf_no_ol_gwas_CEU  = $out_dir_gwas   .  '_temp_nool_LD_PLINK_gwasCEU.vcf';
#my $temp_vcf_no_ol_gwas_YRI  = $out_dir_gwas   .  '_temp_nool_LD_PLINK_gwasYRI.vcf';

#temp: all gwas SNPs not falling in LD blocks, enhanced to look like the output of bedtools intersect -wo
my $temp_bed_no_ol_gwas_CEU  = $out_dir_gwas   .  '_temp_nool_LD_PLINK_gwasCEU.bed';
#my $temp_bed_no_ol_gwas_YRI  = $out_dir_gwas   .  '_temp_nool_LD_PLINK_gwasYRI.bed';

my $temp_bed_ol_gwas_CEU     = $out_dir_gwas   . '_temp_ol_LD_PLINK_gwasCEU.bed';
#my $temp_bed_ol_gwas_YRI     = $out_dir_gwas   . '_temp_ol_LD_PLINK_gwasYRI.bed';

my $out_bed_gwas_CEU     = $out_dir_gwas   .  'LD_PLINK_CEU_GWAScat.maxp_5E-8.bed.gz';
#my $out_bed_gwas_YRI     = $out_dir_gwas   .  'LD_PLINK_YRI_GWAScat.bed.gz';

#output columns of bedtools -wo will be
#CHROM	START	STOP	SNPNAMES	KB_BLOCK	GWAS_CHR	GWAS_POS	GWAS_RSID	[OTHER VCF]	BEDTOOLS_INTERSECT
#so for all those GWAS snps which do not intersect any LD blocks, I will have to append their vcf data to a
#CHROM	START	STOP	SNPNAME=1	KB_BLOCK=-

#temp bed of LD blocks overlapping GWAS snps
system "$BEDTOOLS intersect -wo -a $INPUT_LD_CEU -b $INPUT_GWAS  > $temp_bed_ol_gwas_CEU";
#system "$BEDTOOLS intersect -wo -a $INPUT_LD_YRI -b $INPUT_GWAS  > $temp_bed_ol_gwas_YRI";

#temp vcf of GWAS snps NOT overlapping any LD block
#the format of the output is
#1       7281492 rs766156        -       -       -       Maternal transmission distortion        Gender;Male;Female;Reproductive(p=1.606E-22)
system "$BEDTOOLS intersect -v -a $INPUT_GWAS  -b $INPUT_LD_CEU > $temp_vcf_no_ol_gwas_CEU";
#system "$BEDTOOLS intersect -v -a $INPUT_GWAS  -b $INPUT_LD_YRI > $temp_vcf_no_ol_gwas_YRI";

#prepare full row [LD_BED] + [GWAS_VCF]
#[CHROM	START	STOP	SNPNAMES	KB_BLOCK]	+ [GWAS_CHR	GWAS_POS	GWAS_RSID	[OTHER VCF]	BEDTOOLS_INTERSECT]
#open each temp file, update adding the following columns:
#CHROM	START	STOP	SNPNAME=(1)rsid	KB_BLOCK=0

my $KB_LENGHT='-';

#GWAS - CEU
my %no_ol_gwas_CEU;
open (my $instream,     q{<}, $temp_vcf_no_ol_gwas_CEU) or die("Unable to open $temp_vcf_no_ol_gwas_CEU : $!");
while(<$instream>){
	chomp;
	my @fields = split("\t", $_);
	my $new_bed_line = $fields[0] . "\t" . ( $fields[1]- 1 ) . "\t" . $fields[1] . "\t(1)" . $fields[2] . "\t" .  '-';
	my $new_line = $new_bed_line . "\t" . $_ . "\t" . '0';
	$no_ol_gwas_CEU{$new_line} = 1;
}
close $instream;
open (my $outstream,  q{>}, $temp_bed_no_ol_gwas_CEU) or die("Unable to open $temp_bed_no_ol_gwas_CEU : $!");
foreach my $item (keys %no_ol_gwas_CEU){ print $outstream $item, "\n"; }
close $outstream;

##GWAS - YRI
#my %no_ol_gwas_YRI;
#tie *FILE,   'IO::Zlib', $temp_vcf_no_ol_gwas_YRI, "rb";
##open ($instream,     q{<}, $temp_vcf_no_ol_gwas_YRI) or die("Unable to open $temp_vcf_no_ol_gwas_YRI : $!");
#while(<FILE>){
#	chomp;
#	my @fields = split("\t", $_);
#	my $new_bed_line = $fields[0] . "\t" . ( $fields[1]- 1 ) . "\t" . $fields[1] . "\t(1)" . $fields[2] . "\t" .  '-';
#	my $new_line = $new_bed_line . "\t" . $_ . "\t" . '0';
#	$no_ol_gwas_YRI{$new_line} = 1;
#}
#close FILE;
#open ($outstream,  q{>}, $temp_bed_no_ol_gwas_YRI) or die("Unable to open $temp_bed_no_ol_gwas_YRI : $!");
#foreach my $item (keys %no_ol_gwas_YRI){ print $outstream $item, "\n"; }
#close $outstream;


#cat files, sort, gzip
system "cat $temp_bed_ol_gwas_CEU  $temp_bed_no_ol_gwas_CEU  | sort -k1,1V -k2,2g | gzip -c > $out_bed_gwas_CEU";
#system "cat $temp_bed_ol_gwas_YRI  $temp_bed_no_ol_gwas_YRI  | sort -k1,1V -k2,2g | gzip -c > $out_bed_gwas_YRI";


unlink $temp_vcf_no_ol_gwas_CEU;
#unlink $temp_vcf_no_ol_gwas_YRI;
unlink $temp_bed_no_ol_gwas_CEU;
#unlink $temp_bed_no_ol_gwas_YRI;
unlink $temp_bed_ol_gwas_CEU;
#unlink $temp_bed_ol_gwas_YRI;
print STDERR "FINISHED.\n";