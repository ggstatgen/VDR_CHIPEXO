#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#convert the list of candidate VDR-QTL snps to vcf, so that you can intersect it with other stuff using bedtools

#The MINIMUM vcf understood by bedtools contains the following fields
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

#sample SNPTEST output

#chr peak_id max_peak_rd alternate_ids rsid chromosome position alleleA alleleB index average_maximum_posterior_call info cohort_1_AA cohort_1_AB cohort_1_BB cohort_1_NULL all_AA all_AB all_BB all_NULL all_total all_maf missing_data_proportion bayesian_add_log10_bf bayesian_add_beta_1 bayesian_add_se_1 comment
#chr1 2 21.722534503026 --- rs114983708 NA 714019 A G 2 0.956813 0.831906 20.719 6.25601 0.025026 0 20.719 6.25601 0.025026 0 27 0.116779 0 -0.0557087 0.164558 0.306999 NA

#chrom - field 0
#pos - field 6
#id - field 2
#ref -
#alt -
#QUAL - 
#FILTER -
#INFO - ASSAY, BF, peakID, MAX_PEAK_RD,A,B,AA,AB,BB,BETA,SE

my $INFO_ASSAY = 'VDR-QTL-SNPTEST';
my $INFO_BF_def = 'BF='; #field 
my $INFO_PEAKID_def = 'PEAK_ID='; #field 4
my $INFO_MAX_PEAK_RD_def = 'MAX_PEAK_RD='; #field 5
my $INFO_A_def = 'A=';
my $INFO_B_def = 'B=';
my $INFO_AA_def = 'AA=';
my $INFO_AB_def = 'AB=';
my $INFO_BB_def = 'BB=';
my $INFO_BETA_def = 'BETA=';
my $INFO_SE_def = 'SE=';

my $infile;
my $MIN_BAYES_FACTOR;
GetOptions(
        'i=s'         =>\$infile,
        'minbp=i'     =>\$MIN_BAYES_FACTOR
);

$infile = '/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_diffbind/overlap_3/bayesian_cov_sex_ethn_0.2/snptest_vdr_o3_matrix_p_global_out.txt';

if(!$infile){
     print "USAGE: do_association_SNPTEST_out_to_vcf.pl -i=<INFILE> -minbp=<MIN_BAYES_FACTOR>\n";
     print "<INFILE> SNPTEST output file\n";
     print "(optional)<MIN_BAYES_FACTOR> threshold on BF value (default: 1.0)\n";
     exit 1;
}
if(!$MIN_BAYES_FACTOR){
	$MIN_BAYES_FACTOR = 1;
}

#print header
print "##fileformat=VCFv4.1\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if ($_ =~ /bayesian_add_beta_1/);
	my @fields = split(' ', $_);
	
	next if($fields[23] <  $MIN_BAYES_FACTOR);
	
	my $chr = $fields[0];
	$chr =~ s/chr(.+)/$1/;
	my $pos =  $fields[6];
	my $rsID = $fields[4];
	
	
	#INFO bits of data
	my $INFO_BF = $INFO_BF_def . $fields[23];
	my $INFO_PEAKID =  $INFO_PEAKID_def . $fields[1];
	my $INFO_MAX_PEAK_RD = $INFO_MAX_PEAK_RD_def . $fields[2];
	my $INFO_A    = $INFO_A_def  . $fields[7];
	my $INFO_B    = $INFO_B_def  . $fields[8];
	my $INFO_AA   = $INFO_AA_def . $fields[16];
	my $INFO_AB   = $INFO_AB_def . $fields[17];
	my $INFO_BB   = $INFO_BB_def . $fields[18];
	my $INFO_BETA = $INFO_BETA_def . $fields[24];
	my $INFO_SE   = $INFO_SE_def   . $fields[25];

	my @INFO_ARRAY;
	push @INFO_ARRAY, $INFO_ASSAY, $INFO_BF, $INFO_PEAKID, $INFO_MAX_PEAK_RD, $INFO_A, $INFO_B, $INFO_AA, $INFO_AB, $INFO_BB, $INFO_BETA, $INFO_SE;
	my $INFO = join(";", @INFO_ARRAY); 

	#build data line
	#currently REF and ALT are set to '-'
	#QUAL will contain the PICS score (field 5)
	#FILTER will contain the DISEASE (field 0)
	my $REF = '-';
	my $ALT = '-';
	my $QUAL = '-';
	my $FILTER = '-';
	
	my $data_line = $chr  . "\t" .
					$pos  . "\t" . 
					$rsID . "\t" . 
					$REF  . "\t" .  
	                $ALT  . "\t" .
	                $QUAL . "\t" . 
	                $FILTER   . "\t" .
	                $INFO;
	print $data_line, "\n";	                          
}
close $instream;