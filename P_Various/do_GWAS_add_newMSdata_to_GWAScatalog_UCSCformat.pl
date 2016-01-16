#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#I want to add the beecham data to the UCSC format catalog

#INFO
#There is this Nature Genetics Paper [2013]    (PMED ID 24076602)
#http://www.nature.com/ng/journal/v45/n11/full/ng.2770.html?WT.ec_id=NG-201311
#which has two tables with 48+48 new MS-associated SNPs
#I cut+pasted the tables into a file - now I want to create a "spoof" new GWAS catalog to include these SNPs
#This program receives as input a GWAScatalog.txt file and a tsv file with the new snps
#It outputs a new GWAScatalog_plus_newMS.txt which you then pass to do_GWAS_catalog_to_bed.pl

#GWAS CATALOG FIELDs (UCSC format)
#0 592	
#1 chr1
#2 1005805	
#3 1005806
#4 rs3934834
#5 19851299
#6 Johansson A
#7 2009-10-22
#8 Obesity (Silver Spring)
#9 Linkage and genome-wide association analysis of obesity-related phenotypes: association of weight with the MGAT1 gene.
#10 Body mass index
#11 1,079 South Tyrolian individuals, 790 Dutch founder individuals, 2,060 European ancestry individuals
#12 NA
#13 1p36.33
#14 NR
#15 rs3934834-G
#16 0.80
#17 6E-7
#18 (females + males)
#19 .11
#20 [NR] kg increase
#21 Illumina [318,237]
#22 N


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

#CONSTANTS (fields in the GWAS catalog that are not in the new data)
my $F00_catalog_id = 'NA'; #not sure what this is
my $F01_CHR_ID;   #AVAILABLE
my $F02_START; #$F03_POS - 1
my $F03_POS; #AVAILABLE
my $F04_RSID; #AVAILABLE
my $F05_PUBMEDID = '24076602';
my $F06_FIRSTAUTHOR = 'Beecham AH';
my $F07_DATE = '12/11/2013'; #some date
my $F08_JOURNAL = 'Nat Genet';
my $F09_STUDY = 'Analysis of immune-related loci identifies 48 new susceptibility variants for multiple sclerosis';
my $F10_DISEASE = 'Multiple Sclerosis';
my $F11_SAMPLESIZE = '14,498 subjects and 24,091 controls (european)';
my $F12_REPSIZE = '+14,802 subjects from former GWAS and +26,703 controls (european)';
my $F13 = '';
my $F14_GENE; #available (12) 
my $F15_SNP_RA; #AVAILABLE (0-3)
my $F16_RAF; #AVAILABLE (7?)
my $F17_PVAL; #AVAILABLE
my $F18 = ''; #not sure what this contains; in the GWAS catalog it is listed alongside pval
my $F19_ODDSRATIO; #AVAILABLE (11)
my $F20_CI = 'check_suppl_mat'; #available, but not in the table: it's in a supplementary .xls file. If you find overlaps, add ci from the table
my $F21_platform = 'Illumina ImmunoChip custom [196,524]';
my $F22 = 'N';

my $infile_gwas;
my $infile_ms;
GetOptions(
        'i1=s'      =>\$infile_gwas,
        'i2=s'      =>\$infile_ms
);
if(!$infile_gwas){
     print "USAGE: do_add_newMSdata_to_GWAScatalog_UCSCformat.pl -i1=<INFILE1> -i2=<INFILE2>\n";
     print "<INFILE1> UCSC GWAScatalog .gz file\n";
     print "<INFILE2> MS data tab separated file\n";
     exit 1;
}
if(!$infile_ms){
     print "USAGE: do_add_newMSdata_to_GWAScatalog.pl -i1=<INFILE1> -i2=<INFILE2>\n";
     print "<INFILE1> GWAScatalog txt file\n";
     print "<INFILE2> MS data tab separated file\n";
     exit 1;
}

#open new data and create a new line in the GWAS catalog format; store in hash, then dump all together to stdout
open (my $instream,  q{<}, $infile_ms) or die("Unable to open $infile_ms : $!");
my %gwas_like_line;
while(<$instream>){
	chomp;
	next if($_ eq '');
	($F04_RSID, $F01_CHR_ID, $F03_POS, $F15_SNP_RA, $F16_RAF, $F17_PVAL,$F19_ODDSRATIO,$F14_GENE) = (split /\t/)[0,1,2,3,7,10,11,12];
	
	#field 15 is snp-allele therefore
	$F15_SNP_RA = $F04_RSID . '-' . $F15_SNP_RA;
	$F01_CHR_ID = 'chr' . $F01_CHR_ID; #hg19 output
	
	#the pvalue needs to be reformatted
	#'4.4  1011' > '4.4E-11
	#all these pvalues are < 5E10-8
	#the UCSC format is 5E-6 perl understands it
	#replace '  10' with 'E-' ?
	$F17_PVAL =~ s/\s\s10/E-/;	
	if($F14_GENE eq ''){ $F14_GENE = '-'; }
	
	
	my $line = $F00_catalog_id  . "\t" .
			   $F01_CHR_ID      . "\t" .
			   ($F03_POS-1)     . "\t" .
			   $F03_POS         . "\t" .
			   $F04_RSID        . "\t" .
			   $F05_PUBMEDID    . "\t" .
			   $F06_FIRSTAUTHOR . "\t" . 
			   $F07_DATE        . "\t" .
			   $F08_JOURNAL     . "\t" .
			   $F09_STUDY       . "\t" .
			   $F10_DISEASE     . "\t" .
			   $F11_SAMPLESIZE  . "\t" .
			   $F12_REPSIZE     . "\t" .
			   $F13             . "\t" .
			   $F14_GENE        . "\t" .
			   $F15_SNP_RA      . "\t" .
			   $F16_RAF         . "\t" .
			   $F17_PVAL        . "\t" .
			   $F18             . "\t" .
			   $F19_ODDSRATIO   . "\t" .
			   $F20_CI          . "\t" .
			   $F21_platform    . "\t" .
			   $F22             . "\n";
			   
	$gwas_like_line{$line} = 1;
}
close $instream;

tie *FILE,   'IO::Zlib', $infile_gwas, "rb";
while (<FILE>){
	next if $_ eq '';
	print $_;
}
close FILE;

#now dump the new "gwas catalog extension"
foreach my $item (keys %gwas_like_line){ print $item; }

