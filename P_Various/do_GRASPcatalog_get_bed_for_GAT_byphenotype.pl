#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#2015/4/7
#modified to run with GRASP v 2 to give Ram results with r^2 data

#I have intersected my SNPs from the VDR-QTL or alleleseq analysis with GRASP and found some intersections.
#I want to run a GAT analysis to see if my VDR QTL or whatever snps are intersecting WITH THE SUBSET OF PHENOTYPES OBSERVED more frequently than it would be expected by chance

#notice this does not slop the output, it just returns the snp's coordinate. IT SHOULD overlap for the traits specified when you use this in GAT otherwise you have made an error with the coordinates

#SAMPLE INPUT
#NHLBIkey	HUPfield	LastCurationDate	CreationDate	SNPid(dbSNP134)	chr(hg19)	pos(hg19)	PMID	SNPid(in paper)	LocationWithinPaper	Pvalue	Phenotype	PaperPhenotypeDescription	PaperPhenotypeCategories	DatePub	InNHGRIcat(as of 3/31/12)	Journal	Title	IncludesMale/Female Only Analyses	Exclusively Male/Female	Initial Sample Description	Replication Sample Description	Platform [SNPs passing QC]	GWASancestryDescription	TotalSamples(discovery+replication)	TotalDiscoverySamples	European Discovery	African Discovery	East Asian Discovery	Indian/South Asian Discovery	Hispanic Discovery	Native Discovery	Micronesian Discovery	Arab/ME Discovery	Mixed Discovery	Unspecified Discovery	Filipino Discovery	Indonesian Discovery	Total replication samples	European Replication	African Replication	East Asian Replication	Indian/South Asian Replication	Hispanic Replication	Native Replication	Micronesian Replication	Arab/ME Replication	Mixed Replication	Unspecified Replication	Filipino Replication	Indonesian Replication	InGene	NearestGene	InLincRNA	InMiRNA	InMiRNABS	dbSNPfxn	dbSNPMAF	dbSNPalleles/het/se	dbSNPvalidation	dbSNPClinStatus	ORegAnno	ConservPredTFBS	HumanEnhancer	RNAedit	PolyPhen2	SIFT	LS-SNP	UniProt	EqtlMethMetabStudy

#You want 
#chr(hg19)
#pos-1 (hg19)
#pos (hg19)
#Phenotype

#inputs
#comma-separated list of phenotypes (you will have gotten these from the GRASP catalog so use the same name)


my $phenotype_file;
my $INFILE_CATALOG = "/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/GRASP2final.gz";
GetOptions(
        'ph=s'      =>\$phenotype_file,
);
if(!$phenotype_file){
     print "USAGE: do_GRASPcatalog_get_bed_for_GAT_byphenotype.pl -ph=<FILE>\n";
     print "<FILE> text file with one phenotype per row, including only the phenotypes you want included in the output .bed (phenotypes should be from GRASP V2)\n";
     print "NEEDS IMPROVEMENT! IT WILL ONLY FIND EXACT PHENOTYPE NAMES\n";
     exit 1;
}
#get list
my %phenotypes;
open (my $instream,  q{<}, $phenotype_file) or die("Unable to open $phenotype_file : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	$phenotypes{lc($_)} = 1;
}
close $instream;


tie *FILE,   'IO::Zlib', $INFILE_CATALOG, "rb";
my %unique_entries;
#print "Collecting unique entries for the phenotypes specified..\n";
#for me a unique entry is a variant identified by chr - pos - phen 
#I don't care about pvalues or other id: if there is phenotype I'm interested in at chr - pos, I store it
while (<FILE>)	{ 
	chomp;
	my $data = $_;
	
	next if($_ eq '');
	next if($_ =~ /^NHLBIkey/); #header
	
	my ($key, $snp_id, $chr, $pos, $paper_id, $pval, $phenotype, $paper_phenotype_description) = (split /\t/)[0,4,5,6,8,10,11,12];
	
	#8 aprile 2015
	#final doc has "height" in 11 and "Primary tooth eruption" in 12
	#double check we are referring to this height?
	#TODO TEMP
	if($phenotype =~ /Height/i){
		next unless($paper_phenotype_description =~ /Primary tooth eruption/i);
	}

	next unless($phenotypes{lc($phenotype)});
	
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
	#there are strange items in some chromosome numbers so clean up
	next unless( ($chr =~ /^\d+$/) || ($chr =~ /^[XYM]$/) );
	$chr = 'chr' . $chr;
	
	my $line = $chr . "\t" . ($pos - 1) . "\t" . $pos . "\t" . $phenotype;
	$unique_entries{$line} = 1;
}
close FILE;

foreach my $item (sort keys %unique_entries){
	print $item, "\n";
}