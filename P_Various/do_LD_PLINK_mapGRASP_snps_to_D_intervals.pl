#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#Given the 1000g LD blocks for CEU and YRI obtained with Plink  and stored here
#/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_CEU_LD_b37.bed.gz
#/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_YRI_LD_b37.bed.gz
#/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_CEU_LD_b37_noindels_MAF0.01.bed.gz

#I want to create a bed where I replace, for each GWAS snp reaching genome-wide significance, its LD block.
#crucially, if I find a GWAS snp does not intersect any LD block, the output bed will feature a 1bp interval with that GWAS snp's coordinates

#This will be similar to the DistilLD output I believe.
#many disease snps will correspond to the same LD interval.

#use bedtools as follows to replace the gwas snps INTERSECTING blocks with their block:
#bedtools intersect -wo -a /net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_CEU_LD_b37.bed.gz -b GRASP2final.vcf.gz | gzip -c > LD_PLINK_CEU_GRASP2final.vcf.gz
#bedtools intersect -wo -a /net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_YRI_LD_b37.bed.gz -b GRASP2final.vcf.gz | gzip -c > LD_PLINK_YRI_GRASP2final.vcf.gz

#inputs:
#1 vcf.gz of the gwas catalog obtained with do_GWAScatalog_to_vcf.pl
#2 vcf.gz of the GRASP catalog obtained with do_GRASPcatalog_to_vcf.pl
#3 LD blocks CEU
#4 LD blocks YRI
my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";

my $INPUT_LD;
my $INPUT_GRASP;
GetOptions(
        'i_ld=s'         =>\$INPUT_LD,
        'i_gwas=s'       =>\$INPUT_GRASP
);
if(!$INPUT_LD){
     print "USAGE: do_LD_PLINK_mapGRASP_snps_to_D_intervals.pl -i_ld=<INPUT_LD_FILE> -i_gwas=<INPUT_GRASP_VCF>\n";
     print "<INPUT_LD_FILE> gzipped bed containing D' ld blocks, b37 format\n";
     print "<INPUT_GRASP> gzipped vcf file containing genome wide significant variants from the GRASP catalog\n";
     exit 1;
}
if(!$INPUT_GRASP){
     print "USAGE: do_LD_PLINK_mapGRASP_snps_to_D_intervals.pl -i_ld=<INPUT_LD_FILE> -i_gwas=<INPUT_GRASP_VCF>\n";
     print "<INPUT_LD_FILE> gzipped bed containing D' ld blocks, b37 format\n";
     print "<INPUT_GRASP> gzipped vcf file containing genome wide significant variants from the GRASP catalog\n";
     exit 1;
}
#DISEASE INFO
#for both, only snps with p < 5e-8; only traits with 5+ snps at least
#my $INPUT_GRASP = "/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/GRASP2final.vcf.gz";


#LD INTERVALS
#my $INPUT_LD_CEU = "/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_CEU_LD_b37.bed.gz";
#my $INPUT_LD_YRI = "/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_YRI_LD_b37.bed.gz";

#derive ethnicity
my $ethnicity;
if($INPUT_LD =~ /(ceu)/i){
	$ethnicity = $1;
}elsif($INPUT_LD =~ /(yri)/i){
	$ethnicity = $1;
}else{
	print "Error: impossible to determine ethnicity of LD population used from filename: $INPUT_LD.\nAborting..\n";
	exit -1;
}
#output directory will be where the disease data is 
my ($filename_grasp, $out_dir_grasp) = fileparse($INPUT_GRASP);
print STDERR "Saving all output under $out_dir_grasp..\n";

my $temp_bed_ol    = $out_dir_grasp . '_temp_ol_LD_PLINK_grasp_'   . $ethnicity . '.bed';
my $temp_vcf_no_ol = $out_dir_grasp . '_temp_nool_LD_PLINK_grasp_' . $ethnicity . '.vcf';

my $temp_bed_no_ol = $out_dir_grasp . '_temp_nool_LD_PLINK_grasp_' . $ethnicity . '.bed';
my $out_bed        = $out_dir_grasp .  'LD_PLINK_GRASP_new'           . $ethnicity . '.bed.gz';

system "$BEDTOOLS intersect -wo -a $INPUT_LD -b $INPUT_GRASP > $temp_bed_ol";
system "$BEDTOOLS intersect -v -a $INPUT_GRASP -b $INPUT_LD  > $temp_vcf_no_ol";
#GRASP - CEU
my %no_ol_grasp;
tie *FILE,   'IO::Zlib', $temp_vcf_no_ol, "rb";
while(<FILE>){
	chomp;
	my @fields = split("\t", $_);
	my $new_bed_line = $fields[0] . "\t" . ( $fields[1]- 1 ) . "\t" . $fields[1] . "\t(1)" . $fields[2] . "\t" .  '-';
	my $new_line = $new_bed_line . "\t" . $_ . "\t" . '0';
	$no_ol_grasp{$new_line} = 1;
}
close FILE;
open (my $outstream,  q{>}, $temp_bed_no_ol) or die("Unable to open $temp_bed_no_ol : $!");
foreach my $item (keys %no_ol_grasp){ print $outstream $item, "\n"; }
close $outstream;

system "cat $temp_bed_ol $temp_bed_no_ol| sort -k1,1V -k2,2g | gzip -c > $out_bed";

unlink $temp_vcf_no_ol;
unlink $temp_bed_no_ol;
unlink $temp_bed_ol;
print STDERR "FINISHED.\n";