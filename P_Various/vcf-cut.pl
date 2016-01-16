#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#16/10/2014
#purpose:
#1)SNAP and haploreg are using 1000g pilot data. Can I find more SNPs in LD by manually getting the LD blocks?
#2)I want to obtain LD blocks to do a GAT analysis using the ASB snps and see if they're enriched in particular GWAS intervals
#I found this thread:
#https://www.biostars.org/p/2909/#75824
#and this post, proposing a Make pipeline
#Here is a Makefile.in which I used to generate LD data from 1kG phase1 data. It can do everything from downloading the vcf files to calculation LD using Intersnp or Plink. Some perl scripts and binaries are missing, but that should serve as an example at least. One can use make -j 8 to run the processes in parallel, but make download should have run first.
#However their file vcf-cut.pl is missing. This should try to replace it

#USAGE FOR THE ORIGINAL (MISSING) vcf-cut.pl:
#vcf-cut.pl -c panel genotypes.vcf > genotypes.subset.vcf
#INPUTS:
#1) panel: a text file extracted by the pipeline from the 1000g panel file:
#ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel
#panel contains the subset of samples for a population, eg EUR:
#NA10847
#NA10851
#NA11829
#NA11830
#NA11831
#NA11843
#NA11892
#NA11893
#NA11894
#..
#2) genotypes.vcf.gz vcf/genotypes file from 1000g, eg:
#ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr7.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz

#This script should call the vcf-tools  vcf-subset tool to slice this .vcf.gz to only keep the samples IN PANEL
#so similar to what I do in /net/isi-backup/giuseppe/scripts/do_vcf_slice_by_sample.sh
my $VCFTOOLS = "/net/isi-scratch/giuseppe/tools/vcftools_0.1.12b/bin/vcf-subset";
#sample usage
#${PCODE}/vcf-subset -e -c ${ALLELESEQ_1KG} ${INPUT_VCF} | bgzip -c > ${OUT}/${BASENAME}_ALLELESEQ_1KG.vcf.gz

my $input_vcf;
my $input_panel;
GetOptions(
        'vcf=s'      =>\$input_vcf,
        'panel=s'    =>\$input_panel
);
if(!$input_vcf){
     print "USAGE: vcf-cut.pl -vcf=<VCF_IN> -panel=<PANEL>\n";
     print "<VCF_IN> vcf.gz file from 1000g with genotypes\n";
     print "<PANEL> panel file with subset of samples to keep (eg all CEU, all YRI, etc)\n";
     exit 1;
}
if(!$input_panel){
     print "USAGE: vcf-cut.pl -vcf=<VCF_IN> -panel=<PANEL>\n";
     print "<VCF_IN> vcf.gz file from 1000g with genotypes\n";
     print "<PANEL> panel file with subset of samples to keep (eg all CEU, all YRI, etc)\n";
     exit 1;
}


my ($filename, $directory) = fileparse($input_vcf);
$filename =~ s/(.*)\..*/$1/;
my $outfile = $filename . '.subset.vcf.gz';
#read the panel of samples and turn into comma separated list for input to vcf tools
my @samples;
open (my $instream,  q{<}, $input_panel) or die("Unable to open $input_panel : $!");
while(<$instream>){
	chomp;
	push(@samples, $_);
}
close $instream;
my $sample_list = join(",",@samples);
#system "$VCFTOOLS -e -c $sample_list $input_vcf | bgzip -c > $outfile";
system "$VCFTOOLS -e -c $sample_list $input_vcf";