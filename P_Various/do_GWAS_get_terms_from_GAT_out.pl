#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;


#I want to print a subset of lines corresponding to an input list of phenotypes
#I will need these lines to print an R bar plot to compare them
#this get all the lines  from a GAT .tsv output file corresponding to phenotypes in an input file and returns them in a file
#use it via cluster to run it on all .tsv in a directory

my $infile_gat;
my $infile_phenotypes;
GetOptions(
        'i_gat=s'      =>\$infile_gat,
        'i_ph=s'       =>\$infile_phenotypes
);
if(!$infile_gat){
     print "USAGE: do_GWAS_get_terms_from_GAT_out.pl -i_gat=<GAT_TSV> -i_ph=<PH_LIST>\n";
     print "<INFILE_GAT> GAT .tsv output file\n";
     exit 1;
}