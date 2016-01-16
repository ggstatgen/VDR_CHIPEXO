#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#old prefunseq


#This program accepts 2 input lists of GWAS phenotypes and returns those that intersect
#useful to check which phenotypes are simultaneously 1)mapped to ASB snps replicating by position and 2) such that ASB snps are enriched in the disease LD block compared to a background

#How I generated the final list of diseases for the table in the paper:

#to find the final set of terms, intersect:
#list of terms coming out of the [GRASP/GWAScatalog] analysis as enriched at q < 0.1: I ask the question: if my ASB SNP set intersects LD blocks strongly associated #with DISEASE, for which diseases is this overlap higher than expected by chance?
#INTERSECT THIS WITH:
#set of ASB snps which REPLICATE at the same position in at LEAST 2 samples
#THEN once you have a set of diseases intervals enriched with ASB snps which are supported by at least 2 samples, tag them as follows:
#do these ASB snps happen in VDR peaks?
#do these ASB snps happen in VDR motifs?
#do these enrichment happen for CEU/YRI LD blocks?
#do these ASB snps replicate in more than 2 samples?
#do these diseases come out from GRASP or GWAS or both?


my $input_1;
my $input_2;
my $catalog_file;
GetOptions(
        'i1=s'      =>\$input_1,
        'i2=s'      =>\$input_2,
);
if(!$input_1){
	print "USAGE: do_LD_plink_04_intersect_term_lists.pl -i1=<PH_LIST1> -i2=<PH_LIST2>\n";
	print "<PH_LIST1> and <PH_LIST2> are two txt files with one phenotype per line, both from GRASP or from GWAScat\n";
     exit 1;
}
if(!$input_2){
	print "USAGE: do_LD_plink_04_intersect_term_lists.pl -i1=<PH_LIST1> -i2=<PH_LIST2>\n";
	print "<PH_LIST1> and <PH_LIST2> are two txt files with one phenotype per line, both from GRASP or from GWAScat\n";
     exit 1;
}

my %info;
open (my $instream,  q{<}, $input_1) or die("Unable to open $input_1 : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	$info{lc($_)} .= "A";
}
close $instream;
open ($instream,  q{<}, $input_2) or die("Unable to open $input_2 : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	$info{lc($_)} .= "B";
}
close $instream;


foreach my $item (sort keys %info) {
	if($info{$item} eq 'AB'){
		#print "$item is in: $info{$item}\n";
		print "$item\n";
	}
}
