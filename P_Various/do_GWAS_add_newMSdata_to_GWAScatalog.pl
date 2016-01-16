#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#INFO
#There is this Nature Genetics Paper [2013]    (PMED ID 24076602)
#http://www.nature.com/ng/journal/v45/n11/full/ng.2770.html?WT.ec_id=NG-201311
#which has two tables with 48+48 new MS-associated SNPs
#I cut+pasted the tables into a file - now I want to create a "spoof" new GWAS catalog to include these SNPs
#This program receives as input a GWAScatalog.txt file and a tsv file with the new snps
#It outputs a new GWAScatalog_plus_newMS.txt which you then pass to do_GWAS_catalog_to_bed.pl

#GWAS CATALOG FIELDs:
#0 Date Added to Catalog
#1 PUBMEDID
#2 First Author
#3 Date
#4 Journal
#5 Link
#6 Study
#7 Disease/Trait
#8 Initial Sample Size
#9 Replication Sample Size
#10 Region
#11 Chr_id
#12 Chr_pos
#13 Reported Gene(s)
#14 Mapped_gene
#15 Upstream_gene_id
#16 Downstream_gene_id
#17 Snp_gene_ids
#18 Upstream_gene_distance
#19 Downstream_gene_distance
#20 Strongest SNP-Risk Allele
#21 SNPs
#22 Merged
#23 Snp_id_current
#24 Context
#25 Intergenic
#26 Risk Allele Frequency
#27 p-Value
#28 Pvalue_mlog
#29 p-Value (text)
#30 OR or beta
#31 95% CI (text)
#32 Platform [SNPs passing QC]
#33 CNV

#MSdata fields: (http://www.nature.com/ng/journal/v45/n11/fig_tab/ng.2770_T1.html)
#0 SNP
#1 Chr
#2 Position
#3 RA
#4 RAF
#5 P
#6 OR	
#7 RAF
#8 P
#9 OR
#10 P
#11 OR
#12 Geneb
#13 Function
#I will use 0,1,2,3,13

#CONSTANTS (fields in the GWAS catalog that are not in the new data)
my $F00_DATEADDED = '12/11/2013'; #american format
my $F01_PUBMEDID = '24076602';
my $F02_FIRSTAUTHOR = 'Beecham AH';
my $F03_DATE = '09/29/2013'; #american format
my $F04_JOURNAL = 'Nat Genet';
my $F05_LINK = 'http://www.ncbi.nlm.nih.gov/pubmed/24076602';
my $F06_STUDY = 'Analysis of immune-related loci identifies 48 new susceptibility variants for multiple sclerosis';
my $F07_DISEASE = 'Multiple Sclerosis';
my $F08_SAMPLESIZE = '14,498 subjects and 24,091 controls';
my $F09_REPSIZE = '+14,802 subjects from former GWAS and +26,703 controls';
my $F10_REGION = ''; #unsure
my $F11_CHR_ID;   #AVAILABLE
my $F12_CHR_POS;  #AVAILABLE (hg19/dbSNP 137)
my $F13_14_15_16_17_18_19_GENEINFO = ''; #I leave this blank for now
my $F20_SNP_RA;   #AVAILABLE
my $F21_SNP;      #AVAILABLE
my $F22_MERGED = ''; #unsure
my $F23_SNP_ID;   #AVAILABLE (rsxxxx, ony xxxx)
my $F24_CONTEXT;  #AVAILABLE ('function')
my $F25_INTERGENIC; #bool based on #24
my $F26_RAF = '';      #AVAILABLE, but I have two. Leave blank for now
my $F27_33  = '';     #pvals are available but not sure I need them now.

my $infile_gwas;
my $infile_ms;
GetOptions(
        'i1=s'      =>\$infile_gwas,
        'i2=s'      =>\$infile_ms
);
if(!$infile_gwas){
     print "USAGE: do_add_newMSdata_to_GWAScatalog.pl -i1=<INFILE1> -i2=<INFILE2>\n";
     print "<INFILE1> GWAScatalog txt file\n";
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
	($F21_SNP, $F11_CHR_ID, $F12_CHR_POS, $F20_SNP_RA, $F24_CONTEXT) = (split /\t/)[0,1,2,3,13];
	
	#field 20 is snp-allele therefore
	$F20_SNP_RA = $F21_SNP . '-' . $F20_SNP_RA;
	#field 23 is rsxxxxx without the rs therefore
	$F23_SNP_ID = $F21_SNP; 
	$F23_SNP_ID =~ s/^rs//;
	#field 25 is boolean if field 24 is "intergenic" therefore
	if(lc($F24_CONTEXT) =~ /intergenic/){ 
		$F25_INTERGENIC = '1';
	}else{
		$F25_INTERGENIC = '0';
	}
	
	#create derivative fields: 
	my $line = $F00_DATEADDED . "\t" .
               $F01_PUBMEDID  . "\t" .
               $F02_FIRSTAUTHOR . "\t" .
               $F03_DATE . "\t" .
               $F04_JOURNAL . "\t" .
               $F05_LINK . "\t" . 
               $F06_STUDY . "\t" . 
               $F07_DISEASE . "\t" . 
               $F08_SAMPLESIZE . "\t" .
               $F09_REPSIZE . "\t" .
               $F10_REGION . "\t" .
               $F11_CHR_ID . "\t" .
               $F12_CHR_POS . "\t" .
               $F13_14_15_16_17_18_19_GENEINFO . "\t" .
               $F13_14_15_16_17_18_19_GENEINFO . "\t" .
               $F13_14_15_16_17_18_19_GENEINFO . "\t" .
               $F13_14_15_16_17_18_19_GENEINFO . "\t" .
               $F13_14_15_16_17_18_19_GENEINFO . "\t" .
               $F13_14_15_16_17_18_19_GENEINFO . "\t" .
               $F13_14_15_16_17_18_19_GENEINFO . "\t" .
               $F20_SNP_RA . "\t" . 
               $F21_SNP . "\t" . 
               $F22_MERGED . "\t" .
               $F23_SNP_ID . "\t" .
               $F24_CONTEXT . "\t" .
               $F25_INTERGENIC . "\t" .
               $F26_RAF . "\t" . 
               $F27_33 . "\t" .
               $F27_33 . "\t" .
               $F27_33 . "\t" .
               $F27_33 . "\t" .
               $F27_33 . "\t" .
               $F27_33 . "\t" .
               $F27_33 . "\n";
	
	$gwas_like_line{$line} = 1;
}
close $instream;

open ($instream,  q{<}, $infile_gwas) or die("Unable to open $infile_gwas : $!");
#dump everything on the stdin
while(<$instream>){ next if $_ eq ''; print $_; }
close $instream;
#now dump the new "gwas catalog extension"
foreach my $item (keys %gwas_like_line){ print $item; }

