#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use IO::Zlib;

#Given the 1000g LD blocks for CEU and YRI obtained with Plink and stored here
#/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_CEU_LD_b37.bed.gz
#/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_YRI_LD_b37.bed.gz

#I CREATED THIS SCRIPT AFTER A CONVERSATION WITH CHRIS N
#HE CLAIMS THAT THE NUMBER OF ASB SNPs in the LD block will influence the GAT result.
#I SHOULD REPLACE ALL SNPs in the same block with one SNP randomly placed in the block
#I want to create a bed where I replace, for each bunch of ASB SNPs to test a 1bp interval to be placed e.g. in the middle of the LD block
#crucially, if I find an ASB snp that  does not intersect any LD block, the output bed will feature a 1bp interval with that GWAS snp's coordinates

#Algorithm

#---------x-x-------------x-----------------x---------  ASB 

#-----__________-------__________--------------------- LD INFO

#---------_----------------_----------------_--------- OUTPUT

#bedtools intersect -u -a <LDINFO> -b <ASB-snps>  --> LD block with intersecting SNPs will appear ONCE. Replace, for each, its midpoint
#then, bedtools -v to find the ASB snps which don't intersect
#bedtools intersect -v -a <ASB-snps> -b <LDINFO> -> snps with no overlap into LD blocks

#inputs:
#1 vcf.gz of the interestinghet ASB snps
#3 LD blocks CEU
#4 LD blocks YRI

#LD INTERVALS
my $INPUT_LD_CEU = "/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_CEU_LD_b37.bed.gz";
my $INPUT_LD_YRI = "/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_YRI_LD_b37.bed.gz";
#ASB SNPS
my $INPUT_ASB = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/ENCODE_FILTERED/OVERLAP/d_VCF/interestingHets_consensus_allsites_o1.vcf.gz";

my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin/bedtools";

#output directory will be where the disease data is 
my ($filename,  $out_dir) = fileparse($INPUT_ASB);
print STDERR "Saving all output under $out_dir..\n";
my $temp_bed_ol_CEU      = $out_dir   . '_temp_ol_LD_PLINK_CEU.bed';
my $temp_bed_ol_YRI      = $out_dir   . '_temp_ol_LD_PLINK_YRI.bed';
my $temp2_bed_ol_CEU     = $out_dir   . '_temp2_ol_LD_PLINK_CEU.bed';
my $temp2_bed_ol_YRI     = $out_dir   . '_temp2_ol_LD_PLINK_YRI.bed';

my $temp_vcf_no_ol_CEU  = $out_dir   . '_temp_nool_LD_PLINK_CEU.vcf';
my $temp_vcf_no_ol_YRI  = $out_dir   . '_temp_nool_LD_PLINK_YRI.bed';
my $temp_bed_no_ol_CEU  = $out_dir   . '_temp_nool_LD_PLINK_CEU.vcf';
my $temp_bed_no_ol_YRI  = $out_dir   . '_temp_nool_LD_PLINK_YRI.bed';

my $out_bed_CEU         = $out_dir   .  'VDR_ASB_tagintervals_o1_LD_PLINK_CEU.bed.gz';
my $out_bed_YRI         = $out_dir   .  'VDR_ASB_tagintervals_o1_LD_PLINK_YRI.bed.gz';

###########################
#1. get LD intervals which overlap at least 1 asb snp
###########################
system "$BEDTOOLS intersect -u -a $INPUT_LD_CEU -b $INPUT_ASB  > $temp_bed_ol_CEU";
system "$BEDTOOLS intersect -u -a $INPUT_LD_YRI -b $INPUT_ASB  > $temp_bed_ol_YRI";

#process those two temp files to replace intervals with a 1bp interval in the middle
my %CEU_overlap;
open (my $instream,     q{<}, $temp_bed_ol_CEU) or die("Unable to open $temp_bed_ol_CEU : $!");
while(<$instream>){
	chomp;
	my @fields = split("\t", $_);
	#replace the interval $fields[1] - $fields[2] with a random nucleotide between the two
	#for example, replace with an interval long 1nt close to the stop
	my $stop = $fields[2];
	my $start = $stop -1;
	my $new_bed_line = $fields[0] . "\t" . $start . "\t" . $stop;
	$CEU_overlap{$new_bed_line} = 1;
}
close $instream;
open (my $outstream,  q{>}, $temp2_bed_ol_CEU) or die("Unable to open $temp2_bed_ol_CEU : $!");
foreach my $item (keys %CEU_overlap){ print $outstream $item, "\n"; }
close $outstream;

my %YRI_overlap;
open ($instream,     q{<}, $temp_bed_ol_YRI) or die("Unable to open $temp_bed_ol_YRI : $!");
while(<$instream>){
	chomp;
	my @fields = split("\t", $_);
	#replace the interval $fields[1] - $fields[2] with a random nucleotide between the two
	#for example, replace with an interval long 1nt close to the stop
	my $stop = $fields[2];
	my $start = $stop -1;
	my $new_bed_line = $fields[0] . "\t" . $start . "\t" . $stop;
	$YRI_overlap{$new_bed_line} = 1;
}
close $instream;
open ($outstream,  q{>}, $temp2_bed_ol_YRI) or die("Unable to open $temp2_bed_ol_YRI : $!");
foreach my $item (keys %YRI_overlap){ print $outstream $item, "\n"; }
close $outstream;


##################
#2. get ASB snps which do not overlap any LD interval block
##################
system "$BEDTOOLS intersect -v -a $INPUT_ASB  -b $INPUT_LD_CEU > $temp_vcf_no_ol_CEU";
system "$BEDTOOLS intersect -v -a $INPUT_ASB  -b $INPUT_LD_YRI > $temp_vcf_no_ol_YRI";

#get bed from these and merge with those above
my %CEU_no_overlap;
open ($instream,     q{<}, $temp_vcf_no_ol_CEU) or die("Unable to open $temp_vcf_no_ol_CEU : $!");
while(<$instream>){
	chomp;
	my @fields = split("\t", $_);
	my $new_bed_line = $fields[0] . "\t" . ( $fields[1]- 1 ) . "\t" . $fields[1];
	$CEU_no_overlap{$new_bed_line} = 1;
}
close $instream;
open ($outstream,  q{>}, $temp_bed_no_ol_CEU) or die("Unable to open $temp_bed_no_ol_CEU : $!");
foreach my $item (keys %CEU_no_overlap){ print $outstream $item, "\n"; }
close $outstream;

my %YRI_no_overlap;
open ($instream,     q{<}, $temp_vcf_no_ol_YRI) or die("Unable to open $temp_vcf_no_ol_YRI : $!");
while(<$instream>){
	chomp;
	my @fields = split("\t", $_);
	my $new_bed_line = $fields[0] . "\t" . ( $fields[1]- 1 ) . "\t" . $fields[1];
	$YRI_no_overlap{$new_bed_line} = 1;
}
close $instream;
open ($outstream,  q{>}, $temp_bed_no_ol_YRI) or die("Unable to open $temp_bed_no_ol_YRI : $!");
foreach my $item (keys %YRI_no_overlap){ print $outstream $item, "\n"; }
close $outstream;


#cat bed files, sort, gzip, save
system "cat $temp2_bed_ol_CEU  $temp_bed_no_ol_CEU  | sort -k1,1V -k2,2g | gzip -c > $out_bed_CEU";
system "cat $temp2_bed_ol_YRI  $temp_bed_no_ol_YRI  | sort -k1,1V -k2,2g | gzip -c > $out_bed_YRI";

unlink $temp2_bed_ol_CEU;
unlink $temp2_bed_ol_YRI;
unlink $temp_bed_no_ol_CEU;
unlink $temp_bed_no_ol_YRI;
unlink $temp_bed_ol_CEU;
unlink $temp_bed_ol_YRI;
unlink $temp_vcf_no_ol_CEU;
unlink $temp_vcf_no_ol_YRI;