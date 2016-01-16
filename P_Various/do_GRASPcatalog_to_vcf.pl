#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#INFO
#This converts the GRASP GWAS/eQTL catalog .txt table found here
#http://apps.nhlbi.nih.gov/Grasp/Overview.aspx
#paper
#http://www.ncbi.nlm.nih.gov/pubmed/24931982
#(the SNP data in the catalog have been mapped to )
#All SNP associations in the Full Download version are mapped to the genome [hg19] build and reference SNP database [dbSNP build 134]
#All GWAS included in GRASP Build 1.0 were first published on or before December 31, 2011

#in a vcf file containing ONLY GENOME WIDE SIGNIFICANT ENTRIES (currently 5 x 10-8) for TRAITS HAVING AT LEAST 5 SNPs
#I need it to intersect it with 1000g PLINK LD blocks 

#TODO can you also save the publication in the vcf somewhere?

#script based on do_GWAScatalog_to_vcf.pl

#SAMPLE INPUT
#NHLBIkey	HUPfield	LastCurationDate	CreationDate	SNPid(dbSNP134)	chr(hg19)	pos(hg19)	PMID	SNPid(in paper)	LocationWithinPaper	Pvalue	Phenotype	PaperPhenotypeDescription	PaperPhenotypeCategories	DatePub	InNHGRIcat(as of 3/31/12)	Journal	Title	IncludesMale/Female Only Analyses	Exclusively Male/Female	Initial Sample Description	Replication Sample Description	Platform [SNPs passing QC]	GWASancestryDescription	TotalSamples(discovery+replication)	TotalDiscoverySamples	European Discovery	African Discovery	East Asian Discovery	Indian/South Asian Discovery	Hispanic Discovery	Native Discovery	Micronesian Discovery	Arab/ME Discovery	Mixed Discovery	Unspecified Discovery	Filipino Discovery	Indonesian Discovery	Total replication samples	European Replication	African Replication	East Asian Replication	Indian/South Asian Replication	Hispanic Replication	Native Replication	Micronesian Replication	Arab/ME Replication	Mixed Replication	Unspecified Replication	Filipino Replication	Indonesian Replication	InGene	NearestGene	InLincRNA	InMiRNA	InMiRNABS	dbSNPfxn	dbSNPMAF	dbSNPalleles/het/se	dbSNPvalidation	dbSNPClinStatus	ORegAnno	ConservPredTFBS	HumanEnhancer	RNAedit	PolyPhen2	SIFT	LS-SNP	UniProt	EqtlMethMetabStudy

#You want 
#SNPid(dbSNP134)
#chr(hg19)
#pos(hg19)
#SNPid(in paper)
#Pvalue 
#Phenotype
#PaperPhenotypeCategories (tells if it's GWAS or eQTL) eg Quantitative trait(s)
#fields 4,5,6,8, 10,11,12

my $PROBABILITY = 0.00000005;
my $MINSNPS = 5;

my $infile;
GetOptions(
        'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_GRASPcatalog_to_vcf.pl -i=<INFILE>\n";
     print "<INFILE> GRASP catalog .gz input file\n";
     print "NOTE: I will only consider SNPs having p < $PROBABILITY and traits with at least $MINSNPS SNPs\n";
     print "CHANGE if needed\n";
     exit 1;
}
#print header
print "##fileformat=VCFv4.1\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

tie *FILE,   'IO::Zlib', $infile, "rb";
#open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
my %unique_entries;
while(<FILE>){
	chomp;
	my $data = $_;
	
	next if($_ eq '');
	next if($_ =~ /^NHLBIkey/); #header
	
	my ($key, $snp_id, $chr, $pos, $pubmed_id, $pval, $phenotype, $ph_category) = (split /\t/)[0,4,5,6,7,10,11,13];
	my $missing = '-'; #REF     ALT     QUAL   
	next if(!$key);
	next if($key eq '');
	next if(!$snp_id);
	next if($snp_id eq '');
	next if(!$chr);
	next if($chr eq '');
	next if(!$pos);
	next if($pos eq '');	
	next if(!$pval);
	next if($pval eq '');
	next if(!$phenotype);
	next if($phenotype eq '');
	
	next unless($pval);
		
	#there are strange items in some chromosome numbers so clean up
	next unless( ($chr =~ /^\d+$/) || ($chr =~ /^[XYM]$/) );
	#filter out QTLs??
	next if($ph_category =~ /Quantitative trait/);
	#filter other stuff
	if($phenotype =~ /^Gene expression of/i){ next;}
	if($phenotype =~ /gene expression in/i){ next;}
	if($phenotype =~ /^Methylation levels/i){ next;}
	if($phenotype =~ /^Differential exon level expression/i){ next;}
	if($phenotype =~ /^Methylation QTL/i){ next;}
	if($phenotype =~ /^PC/){ next;}
	if($phenotype =~ /^SM/){ next;}
	if($phenotype =~ /^Ratio of/){ next;}
	if($phenotype =~ /^Serum metabolite/){ next;}
	if($phenotype =~ /^Serum ratio/){ next;}
	if($phenotype =~ /^Transcript termination/){ next;}
	if($phenotype =~ /^lysoPC/){ next;}
	#next if($phenotype =~ /\(/);
	
	#as of version 2, there are typos in Grave's/Graves' disease and Behcet/ Beh#cet Disease
	$phenotype = "Behcet's disease"   if(lc($phenotype) =~ /^beh(.*)et's disease/);
	$phenotype = "Graves' disease"   if(lc($phenotype) =~ /^grave's disease/);
	$phenotype = "Irritable bowel syndrome"   if(lc($phenotype) =~ /^irritible bowel syndrome/);
	$phenotype = "Psoriasis"   if(lc($phenotype) =~ /^psoriasis (type 1)/);
	$phenotype = "Psoriasis"   if(lc($phenotype) =~ /^psoriasis (cutaneous psoriasis)/);
	
	
	
	#"3.15E?16"
	$pval =~ s/\?/-/;
	next if($pval > $PROBABILITY);
	
	my $filter = $phenotype;
	my $info = $ph_category . '(p=' . $pval . ')';
	$snp_id = 'rs' . $snp_id;
	
	#my $line = $chr . "\t" . $pos . "\t" . $snp_id . "\t" . $missing . "\t" . $missing .  "\t" . $missing . "\t" . $filter . "\t" . $info . "\t". $pubmed_id . "\n";
	#used the above just get the total number of final studies considered
	my $line = $chr . "\t" . $pos . "\t" . $snp_id . "\t" . $missing . "\t" . $missing .  "\t" . $missing . "\t" . $filter . "\t" . $info . "\n";
	$unique_entries{$line} = 1;
}
close FILE;

#2nd pass
#build phenotype->snp-number map based on phenotype
my %phenotype_map;
foreach my $entry (keys %unique_entries){
	my @fields = split (/\t/, $entry);
	if(!$phenotype_map{lc($fields[6])}){
		$phenotype_map{lc($fields[6])} = 1;
	}else{
		$phenotype_map{lc($fields[6])} += 1;		
	}
}

foreach my $item (keys %unique_entries){
	my @fields = split (/\t/, $item);
	next if($phenotype_map{lc($fields[6])} < $MINSNPS); 
	print $item; 
}