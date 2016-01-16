#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#I want to add the beecham data to the GRASP v2 catalog vcf version

#Beecham MSdata fields: (http://www.nature.com/ng/journal/v45/n11/fig_tab/ng.2770_T1.html)
#0 SNP
#1 Chr
#2 Position
#3 RA 
#4 RAF (discovery)
#5 P (discovery)
#6 OR (discovery)	 
#7 RAF (repl)
#8 P (repl)
#9 OR (repl)
#10 P (joint)
#11 OR (joint)
#12 Gene
#13 Function
#I want 0,1,2,3,7(RISK ALLELE FREQ REPL),10,11,12,13

#GRASP catalog vcf
###fileformat=VCFv4.1
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#1       807512  rs10751454      -       -       -       Cardiovascular disease prevalence       Cardiovascular disease (CVD);Myocardial infarction (MI);Mortality;Cancer(p=2.20E-13)

#FILTER: phenotype
#INFO: phenotype category, pvalue

my $infile_grasp_vcf;
my $infile_ms;
GetOptions(
        'i1=s'      =>\$infile_grasp_vcf,
        'i2=s'      =>\$infile_ms
);
if(!$infile_grasp_vcf){
     print "USAGE: do_GRASPcatalog_add_newMSdata.pl -i1=<INFILE1> -i2=<INFILE2>\n";
     print "<INFILE1> GRASP catalog vcf.gz file\n";
     print "<INFILE2> MS data tab separated file\n";
     exit 1;
}
if(!$infile_ms){
     print "USAGE: do_GRASPcatalog_add_newMSdata.pl -i1=<INFILE1> -i2=<INFILE2>\n";
     print "<INFILE1> GRASP catalog vcf.gz file\n";
     print "<INFILE2> MS data tab separated file\n";
     exit 1;
}

#open MS data and create a new line in the vcf format; store in hash, then dump all together to stdout
open (my $instream,  q{<}, $infile_ms) or die("Unable to open $infile_ms : $!");
my %vcf_line;
while(<$instream>){
	chomp;
	next if($_ eq '');
	my ($F04_RSID, $F01_CHR_ID, $F03_POS, $F15_SNP_RA, $F16_RAF, $F17_PVAL,$F19_ODDSRATIO,$F14_GENE) = (split /\t/)[0,1,2,3,7,10,11,12];
	
	#field 15 is snp-allele therefore
	$F15_SNP_RA = $F04_RSID . '-' . $F15_SNP_RA;
	#$F01_CHR_ID = 'chr' . $F01_CHR_ID; #hg19 output
	
	#the pvalue needs to be reformatted
	#'4.4  1011' > '4.4E-11
	#all these pvalues are < 5E10-8
	#the UCSC format is 5E-6 perl understands it
	#replace '  10' with 'E-' ?
	$F17_PVAL =~ s/\s\s10/E-/;	
	if($F14_GENE eq ''){ $F14_GENE = '-'; }
	
#CHROM
#POS
#ID
#REF
#ALT
#QUAL
#FILTER
#INFO

#1       807512  rs10751454      -       -       -       Cardiovascular disease prevalence       Cardiovascular disease (CVD);Myocardial infarction (MI);Mortality;Cancer(p=2.20E-13)
	
	my $info = "Multiple sclerosis, Beecham 2013 (" . $F17_PVAL . ')';
	
	my $line =  $F01_CHR_ID . "\t" . 
				$F03_POS    . "\t" .
				$F04_RSID   . "\t" . 
				"-"         . "\t" .
				"-"         . "\t" .
				"-"         . "\t" .
				'Multiple sclerosis' . "\t" .
				$info . "\n";
				 			
	$vcf_line{$line} = 1;
}
close $instream;

tie *FILE,   'IO::Zlib', $infile_grasp_vcf, "rb";
while (<FILE>){
	next if $_ eq '';
	print $_;
}
close FILE;

#now dump the new "gwas catalog extension"
foreach my $item (keys %vcf_line){ print $item; }
