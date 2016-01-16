#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use IO::Zlib;

#Given the 1000g LD blocks for CEU and YRI obtained with MIGLLD and stored here
#/net/isi-mirror/1000_Genomes/POPULATION_PANELS/LDEXPL_HAPLOTYPE_BLOCKS_r2

#typical structure of the beds:
#chr start end block_id number_of_markers_in_block
#1	61987	61989	BLOCK_0069423	2
#1	86028	86065	BLOCK_0061362	2

#I want to create a bed where I replace, for each GWAS snp reaching genome-wide significance, its LD block.
#crucially, if I find a GWAS snp does not intersect any LD block, the output bed will feature a 1bp interval with that GWAS snp's coordinates

#many disease snps will correspond to the same LD interval.

#inputs:
#2 vcf.gz of the GRASP catalog obtained with do_GRASPcatalog_to_vcf.pl
#3 LD blocks CEU

#LD INTERVALS
#my $INPUT_LD_CEU = "/net/isi-mirror/1000_Genomes/POPULATION_PANELS/LDEXPL_HAPLOTYPE_BLOCKS_r2/MIGLD_CEU_r2_0.8_MAF_0.01_frac_0.9_b37.bed.gz";
my $INPUT_LD_CEU = "/net/isi-mirror/1000_Genomes/POPULATION_PANELS/LDEXPL_HAPLOTYPE_BLOCKS_r2/LD_r0.8_MIG_EUR_noindels_maf0.01_hwe0.001_b37.bed.gz";

#DISEASE INFO
#only snps with p < 5e-8; only traits with 5+ snps at least
#18/2/2015 now using up to date UCSC GWAS cat
#my $INPUT_GWAS = "/net/isi-scratch/giuseppe/indexes/GWAS/UCSC_gwas_catalog/gwasCatalog.b37.maxp_5E-8.minsnp5.vcf.gz";

my $INPUT_GRASP = "/net/isi-scratch/giuseppe/indexes/GWAS_GRASP/GRASP2final_plus_Beecham2013_MS.vcf.gz";
my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";

print STDERR "do_LD_gwassnp2block_MIGLD.pl:\n";
print STDERR "Using the following files:\n";
print STDERR "GRASP catalog vcf: $INPUT_GRASP\n";
print STDERR "LD blocks: $INPUT_LD_CEU\n";
print STDERR "EDIT if required\n";

#output directory will be where the disease data is 
my ($filename_grasp, $out_dir_grasp) = fileparse($INPUT_GRASP);
print STDERR "Saving all output under $out_dir_grasp..\n";

#temp: all gwas snps not falling in LD blocs
my $temp_vcf_no_ol_grasp_CEU = $out_dir_grasp   .  '_temp_nool_LD_MIGLD_graspCEU.vcf';
#my $temp_vcf_no_ol_grasp_YRI = $out_dir_grasp   .  '_temp_nool_LD_PLINK_graspYRI.vcf';

#temp: all gwas SNPs not falling in LD blocks, enhanced to look like the output of bedtools intersect -wo
my $temp_bed_no_ol_grasp_CEU = $out_dir_grasp   .  '_temp_nool_LD_MIGLD_graspCEU.bed';
#my $temp_bed_no_ol_grasp_YRI = $out_dir_grasp   .  '_temp_nool_LD_PLINK_graspYRI.bed';

my $temp_bed_ol_grasp_CEU    = $out_dir_grasp   . '_temp_ol_LD_MIGLD_graspCEU.bed';
#my $temp_bed_ol_grasp_YRI    = $out_dir_grasp   . '_temp_ol_LD_PLINK_graspYRI.bed';

my $out_bed_grasp_CEU    = $out_dir_grasp   .  'LD_r0.8_MIG_EUR_noindels_maf0.01_hwe0.001_GRASP.bed.gz';
#my $out_bed_grasp_YRI    = $out_dir_grasp   .  'LD_PLINK_YRI_GRASP.bed.gz';

#output columns of bedtools -wo will be
#CHROM	START	STOP	BLOCK_ID	MARKERS_IN_BLOCK	GWAS_CHR	GWAS_POS	GWAS_RSID	[OTHER VCF]	BEDTOOLS_INTERSECT
#so for all those GWAS snps which do not intersect any LD blocks, I will have to append their vcf data to a
#CHROM	START	STOP	BLOCK_ID=NA	MARKERS_IN_BLOCK=NA

#temp bed of LD blocks overlapping GWAS snps
system "$BEDTOOLS intersect -wo -a $INPUT_LD_CEU -b $INPUT_GRASP > $temp_bed_ol_grasp_CEU";
#system "$BEDTOOLS intersect -wo -a $INPUT_LD_YRI -b $INPUT_GRASP > $temp_bed_ol_grasp_YRI";

#temp vcf of GWAS snps NOT overlapping any LD block
#the format of the output is
#1       7281492 rs766156        -       -       -       Maternal transmission distortion        Gender;Male;Female;Reproductive(p=1.606E-22)
system "$BEDTOOLS intersect -v -a $INPUT_GRASP -b $INPUT_LD_CEU > $temp_vcf_no_ol_grasp_CEU";
#system "$BEDTOOLS intersect -v -a $INPUT_GRASP -b $INPUT_LD_YRI > $temp_vcf_no_ol_grasp_YRI";

#prepare full row [LD_BED] + [GWAS_VCF]
#[CHROM	START	STOP	BLOCK_ID	MARKERS]	+ [GWAS_CHR	GWAS_POS	GWAS_RSID	[OTHER VCF]	BEDTOOLS_INTERSECT]
#open each temp file, update adding the following columns:
#CHROM	START	STOP	SNPNAME=(1)rsid	KB_BLOCK=0

#GRASP - CEU
my %no_ol_grasp_CEU;
tie *FILE,   'IO::Zlib', $temp_vcf_no_ol_grasp_CEU, "rb";
while(<FILE>){
	chomp;
	my @fields = split("\t", $_);
	my $new_bed_line = $fields[0] . "\t" . ( $fields[1]- 1 ) . "\t" . $fields[1] . "\t(1)" . $fields[2] . "\t" .  '-';
	my $new_line = $new_bed_line . "\t" . $_ . "\t" . '0';
	$no_ol_grasp_CEU{$new_line} = 1;
}
close FILE;
open (my $outstream,  q{>}, $temp_bed_no_ol_grasp_CEU) or die("Unable to open $temp_bed_no_ol_grasp_CEU : $!");
foreach my $item (keys %no_ol_grasp_CEU){ print $outstream $item, "\n"; }
close $outstream;

##GRASP - YRI
#my %no_ol_grasp_YRI;
#tie *FILE,   'IO::Zlib', $temp_vcf_no_ol_grasp_YRI, "rb";
##open ($instream,     q{<}, $temp_vcf_no_ol_grasp_YRI) or die("Unable to open $temp_vcf_no_ol_grasp_YRI : $!");
#while(<FILE>){
#	chomp;
#	my @fields = split("\t", $_);
#	my $new_bed_line = $fields[0] . "\t" . ( $fields[1]- 1 ) . "\t" . $fields[1] . "\t(1)" . $fields[2] . "\t" .  '-';
#	my $new_line = $new_bed_line . "\t" . $_ . "\t" . '0';
#	$no_ol_grasp_YRI{$new_line} = 1;
#}
#close FILE;
#open ($outstream,  q{>}, $temp_bed_no_ol_grasp_YRI) or die("Unable to open $temp_bed_no_ol_grasp_YRI : $!");
#foreach my $item (keys %no_ol_grasp_YRI){ print $outstream $item, "\n"; }
#close $outstream;


#cat files, sort, gzip
#system "cat $temp_bed_ol_gwas_CEU  $temp_bed_no_ol_gwas_CEU  | sort -k1,1V -k2,2g | gzip -c > $out_bed_gwas_CEU";
#system "cat $temp_bed_ol_gwas_YRI  $temp_bed_no_ol_gwas_YRI  | sort -k1,1V -k2,2g | gzip -c > $out_bed_gwas_YRI";
system "cat $temp_bed_ol_grasp_CEU $temp_bed_no_ol_grasp_CEU | sort -k1,1V -k2,2g | gzip -c > $out_bed_grasp_CEU";
#system "cat $temp_bed_ol_grasp_YRI $temp_bed_no_ol_grasp_YRI | sort -k1,1V -k2,2g | gzip -c > $out_bed_grasp_YRI";


unlink $temp_vcf_no_ol_grasp_CEU;
#unlink $temp_vcf_no_ol_grasp_YRI;
#unlink $temp_vcf_no_ol_gwas_CEU;
#unlink $temp_vcf_no_ol_gwas_YRI;
unlink $temp_bed_no_ol_grasp_CEU;
#unlink $temp_bed_no_ol_grasp_YRI;
#unlink $temp_bed_no_ol_gwas_CEU;
#unlink $temp_bed_no_ol_gwas_YRI;
unlink $temp_bed_ol_grasp_CEU;
#unlink $temp_bed_ol_grasp_YRI;
#unlink $temp_bed_ol_gwas_CEU;
#unlink $temp_bed_ol_gwas_YRI;
print STDERR "FINISHED.\n";