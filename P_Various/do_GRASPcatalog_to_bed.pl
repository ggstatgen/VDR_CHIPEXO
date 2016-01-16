#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $BEDTOOLS  = '/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin';
my $CHROM_SIZE_PATH = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes';

#INFO
#This converts the GRASP GWAS/eQTL catalog .txt table found here
#http://apps.nhlbi.nih.gov/Grasp/Overview.aspx
#paper
#http://www.ncbi.nlm.nih.gov/pubmed/24931982
#(the SNP data in the catalog have been mapped to )
#All SNP associations in the Full Download version are mapped to the genome [hg19] build and reference SNP database [dbSNP build 134]
#All GWAS included in GRASP Build 1.0 were first published on or before December 31, 2011

#in a bed file to use for GAT annotation

#script based on do_GWAScatalog_to_bed.pl

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

my $infile;
my $min_snp_number;
my $SLOP;
GetOptions(
        'i=s'      =>\$infile,
        'minsnp=i' =>\$min_snp_number,
        'slop=i'   =>\$SLOP
);

if(!$infile){
     print "USAGE: do_GRASPcatalog_to_bed.pl -i=<INFILE> -minsnp=<MIN_SNP_NUMBER> -slop=<BASEPAIRS>\n";
     print "<INFILE> GRASPcatalog file\n";
     print "<MIN_SNP_NUMBER> minimum number of SNP per phenotype to consider for output\n";
     print "(optional)<BASEPAIRS> number of bp requested upstream (downstream) of each SNP (eg: 150,000)\n";
     exit 1;
}
if(!$min_snp_number){
     print "USAGE: do_GRASPcatalog_to_bed.pl -i=<INFILE> -minsnp=<MIN_SNP_NUMBER> -slop=<BASEPAIRS>\n";
     print "<INFILE> GRASPcatalog file\n";
     print "<MIN_SNP_NUMBER>  minimum number of SNP per phenotype to consider for output\n";
     print "(optional)<BASEPAIRS> number of bp requested upstream (downstream) of each SNP (eg: 150,000)\n";
     exit 1;
}
my ($filename, $directory) = fileparse($infile);
$filename =~ s/(.*)\..*/$1/;
my $outfile_segments = $directory  .  $filename . '_segments.bed';

#1st pass
#collect unique entries (study,disease/trait,chromosome,position,snp id) in hash (in case  there are duplicate SNP-phenotype lines.)
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
my %unique_entries;
print "Collecting unique entries..\n";
while(<$instream>){
	chomp;
	my $data = $_;
	
	next if($_ eq '');
	next if($_ =~ /^NHLBIkey/); #header
	
	my ($key, $snp_id, $chr, $pos, $paper_id, $pval, $phenotype, $ph_category) = (split /\t/)[0,4,5,6,8,10,11,13];
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
	#filter out QTLs??
	next if($ph_category =~ /Quantitative trait/);
	#filter other stuff
	if($phenotype =~ /^gene expression/i){ next;}
	#next if($phenotype =~ /\(/);
	$snp_id = 'rs' . $snp_id;
	$chr = 'chr' . $chr;
	
	my $line = $key . "\t" . $phenotype . "\t" . $chr . "\t" . $pos . "\t" .  $snp_id . "\t" . $pval;
	$unique_entries{$line} = 1;
}
close $instream;

print "Building phenotype->snp-number map..\n";
my %phenotype_map;
foreach my $entry (keys %unique_entries){
	my ($key, $phenotype, $chr, $pos, $snp_id, $pval) = split (/\t/, $entry);
	if(!$phenotype_map{lc($phenotype)}){
		$phenotype_map{lc($phenotype)} = 1;
	}else{
		$phenotype_map{lc($phenotype)} += 1;		
	}
}

my $outfile_phenotypes = $directory . "\/" .  $filename . '_phenotypes.txt';
open (my $temp_stream,  q{>}, $outfile_phenotypes) or die("Unable to open $outfile_phenotypes : $!");
foreach my $item (sort keys %phenotype_map){
	print $temp_stream $item, "\n" if($phenotype_map{$item} >= $min_snp_number);
}

print "Creating per-phenotype bed files and appending to global bed file..\n";
#create one bed for each of the terms above, if the term has > $min_snp_number
my %collective_bed;
foreach my $item (sort keys %phenotype_map){
	next if($phenotype_map{$item} < $min_snp_number);
	
	my $outfile_single_phenotype           = $directory . "\/" .  $filename . '_phenotype_temp.bed';
	my $outfile_single_phenotype_processed = $directory . "\/" .  $filename . '_phenotype_temp_out.bed';
	open (my $temp_stream,  q{>}, $outfile_single_phenotype) or die("Unable to open $outfile_single_phenotype : $!");
	foreach my $entry (keys %unique_entries){
		my ($key, $phenotype, $chr, $pos, $snp_id, $pval) = split (/\t/, $entry);
		#use 0-based start coordinates, which are the proper BED format (pos-1, pos)
		my $line = $chr . "\t" . ($pos-1)  . "\t" . $pos . "\t" . $item . "\n";
		print $temp_stream 	$line;
	}
	close $temp_stream;
	
	#you have a bed with only the snp and its  location. Extend location to +/- 150.000, sort, merge
	if($SLOP){
		system "$BEDTOOLS/bedtools slop -b $SLOP -i $outfile_single_phenotype -g $CHROM_SIZE_PATH |  sort -k1,1V -k2,2g | $BEDTOOLS/bedtools merge -i stdin > $outfile_single_phenotype_processed";		
	}else{
		system "sort -k1,1V -k2,2g $outfile_single_phenotype | $BEDTOOLS/bedtools merge -i stdin > $outfile_single_phenotype_processed";
	}

	open (my $temp_stream_processed,  q{<}, $outfile_single_phenotype_processed) or die("Unable to open $outfile_single_phenotype_processed : $!");
	while(<$temp_stream_processed>){
		chomp;
		next if($_ eq '');
		my $data_line = $_ . "\t" . $item . "\n";
		$collective_bed{$data_line} = 1;
	}
	close $temp_stream_processed;
	unlink $outfile_single_phenotype;
	unlink $outfile_single_phenotype_processed;
}
open (my $segment_stream,  q{>}, $outfile_segments) or die("Unable to open $outfile_segments : $!");
foreach my $item (sort keys %collective_bed){ print $segment_stream $item; }
close $segment_stream;

print "FINISHED\n";