#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;
use List::Util qw(min max sum reduce);

#25/11/2015
#Given the Output.vcf file from funseq, get an .vcf with only the rBVs (the VDR-BVs that reproduce by position in at least 2 samples)

#if you want to remove MHC variants and LD-dependent variants, use do_funseq_VET_01_LDclump.pl

my $input_alleleseq_raw;
my $input_funseq_vcf;

GetOptions(
        'i=s'      =>\$input_funseq_vcf
);
#$INPUT_ALLELESEQ_RAW = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/_RAW/interestingHets_vdrchipexo.txt"; #b37

my $USAGE = "USAGE: do_funseq_get_VDRrBVs.pl -i=<FUNSEQ_OUT>\n". 
            "<FUNSEQ_OUT> .vcf file with the funseq variants to use\n";

if(!$input_funseq_vcf){
	print $USAGE;
    exit 1;
}