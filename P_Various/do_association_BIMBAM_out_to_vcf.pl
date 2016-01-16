#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#convert the list of candidate VDR-QTL snps to vcf, so that you can intersect it with other stuff using bedtools

#The MINIMUM vcf understood by bedtools contains the following fields
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

#sample BIMBAM output
#chr	pos	rsnum	bf	pval	mu	a	d
#7	138803023	rs17160772     	1.851	0.0001	-0.189	0.38	0.015

#chrom - field 0
#pos - field 1
#id - field 2
#ref - 
#alt -
#QUAL -
#FILTER -
#INFO - BF, PVAL,MU,A,D

my $INFO_ASSAY = 'VDR-QTL-BIMBAM';
my $INFO_BF_def = 'BF='; #field 3
my $INFO_PVAL_def = 'PVAL='; #field 4
my $INFO_MU_def = 'MU='; #field 5
my $INFO_A_def = 'A=';
my $INFO_D_def = 'D=';

my $infile;
my $MIN_BAYES_FACTOR;
GetOptions(
        'i=s'         =>\$infile,
        'minbp=i'     =>\$MIN_BAYES_FACTOR
);

$infile = '/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_diffbind/overlap_3/d_BIMBAM/BIMBAM_vdr_o3_matrix_bimbam_global_out.txt';

if(!$infile){
     print "USAGE: do_association_BIMBAM_out_to_vcf.pl -i=<INFILE> -minbp=<MIN_BAYES_FACTOR>\n";
     print "<INFILE> BIMBAM output file\n";
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
	next if ($_ =~ /^chr/);
	
	my $data = $_;
	$data =~ s/ //g;
	my @fields = split("\t", $data);
	
	next if(abs($fields[3]) <  $MIN_BAYES_FACTOR);
	
	my $chr = $fields[0];
	$chr =~ s/chr(.+)/$1/;
	my $pos =  $fields[1];
	my $rsID = $fields[2];
	
	
	#INFO bits of data
	my $INFO_BF = $INFO_BF_def . $fields[3];
	my $INFO_PVAL = $INFO_PVAL_def . $fields[4];
	my $INFO_MU = $INFO_MU_def . $fields[5];
	my $INFO_A = $INFO_A_def . $fields[6];
	my $INFO_D = $INFO_D_def . $fields[7];

	my @INFO_ARRAY;
	push @INFO_ARRAY, $INFO_ASSAY, $INFO_BF, $INFO_PVAL,$INFO_MU,$INFO_A,$INFO_D;
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