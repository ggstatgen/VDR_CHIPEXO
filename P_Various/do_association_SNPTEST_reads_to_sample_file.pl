#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#NOT USED
#see answer on mailing list OXSTATGEN
#Hi Guiseppe,
#No hard threshold that I know of, but SNPTEST is not really designed to work well in that setting at the moment 
#(e.g. currently it'll print out 8000 phenotype summaries in the log file for each scan.)
#So you may be better of using a script to generate per-phenotype sample files and running using those.
#Best wishes,
#Gavin.


#13/5/2014
#This is needed to obtain a sample file for SNPTEST containing the phenotypes (peak sizes)
#sample file format is detailed here:
#http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html#Sample_File_Format_

#INPUT:
#matrix of counts (obtained with do_htseq_get_matrix_from.counts.pl)
#try both normalised and not normalised
#covariate data (sex, ethnicity)

#OUTPUT:
#format follows this convention:
#ID_1 ID_2 missing cov_1 cov_2 cov_3 cov_4 pheno1 bin1
#0 0 0 D D C C P B

#with ID_1, ID_2 and missing compulsory (and corresponding 0 0 0 on second header line) followed by covariates, phenotypes, and binary
#one line per individual

#--------------------------------------
#at the moment I generate the following:
#--------------------------------------
#ID_1 ID_2 missing sex ethnicity peak_1 ... peak_7220
#0 0 0 D D P ... P
#NA06986 NA06986 0.0 male ceu 5.666 ... 12.6
#NA06989 NA06989 0.0 male ceu 4.554 ... 3.45
#NA06997 NA06997 0.0 female yri 2.4345 ... 4.55

#items need to be separated by SPACES, not tabs

#MISSING VALUES: code with NA

my $infile_rcmatrix;
GetOptions(
        'rc=s'          =>\$infile_rcmatrix
);
if(!$infile_rcmatrix){
     print "USAGE: do_ASSOCIATION_bimbam.pl -vcf=<VCF_FILE> -bed=<BED_FILE> -rc=<READ_COUNT>\n";
     print "<VCF_FILE> initial set of variation_genotypes to use(eg hapmap, omni2.5, 1000g)\n";
     print "<BED_FILE> set of peaks to use for the association analysis\n";
     print "<READ_COUNT> matrix of normalised read counts for the peaks in <BED_FILE> + headers. row: peak. column: sample\n";
     exit 1;
}

#eg
$infile_rcmatrix = '/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_htseq_counts/d_all_30/MATRIX.txt';

#############
#1. get the matrix of read counts in a data structure indexed by peak and sample
# also get a hash of candidate sample names (a subset of those in the vcf for which you want to compute regression)
#############
my %peak_to_counts;
my @sample_names;
my %candidate_sample_names; #samples for which we have read counts, could be a subset of those in the vcf
open (my $instream,     q{<}, $infile_rcmatrix) or die("Unable to open $infile_rcmatrix : $!");
while(<$instream>){
	chomp;
	my %sample_name_to_rd;
	if($_ =~ /^\#/){
		@sample_names = split("\t", $_);
		shift @sample_names;
		foreach (@sample_names){ s/\"//g; }
		%candidate_sample_names = map { $_ => 1 } @sample_names;
		next;
	}
	my @counts = split("\t", $_);
	my $this_peak_ID = shift @counts;
	$this_peak_ID =~ s/\"//g;
	@sample_name_to_rd{@sample_names} = @counts;
	
	foreach my $s_name (keys %sample_name_to_rd){
		$peak_to_counts{$this_peak_ID}{$s_name} = $sample_name_to_rd{$s_name};
	}
}
close $instream;




