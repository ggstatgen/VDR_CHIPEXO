#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#INFO 14/11
#This converts the GWAS catalog .txt table found here
#http://www.genome.gov/gwastudies/ (The SNP data in the catalog have been mapped to dbSNP Build 132 and Genome Assembly, GRCh37/hg19. )
#in a minimal vcf to be used with bedtools (I want to intersect it with my peaks)
#The MINIMUM vcf understood by bedtools contains the following fields
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#So create these

#SAMPLE INPUT
#Date Added to Catalog	PUBMEDID	First Author	Date	Journal	Link	Study	Disease/Trait	Initial Sample Size	Replication Sample Size	Region	Chr_id	Chr_pos	Reported Gene(s)	Mapped_gene	Upstream_gene_id	Downstream_gene_id	Snp_gene_ids	Upstream_gene_distance	Downstream_gene_distance	Strongest SNP-Risk Allele	SNPs	Merged	Snp_id_current	Context	Intergenic	Risk Allele Frequency	p-Value	Pvalue_mlog	p-Value (text)	OR or beta	95% CI (text)	Platform [SNPs passing QC]	CNV
#11/13/2013	23696099	Ding K	05/20/2013	G3 (Bethesda)	http://www.ncbi.nlm.nih.gov/pubmed/23696099	Genetic variants that confer resistance to malaria are associated with red blood cell traits in African-Americans: an electronic medical record-based genome-wide association study.	Red blood cell traits	1904 Afican American individuals	411 Afican American individuals	11p15.4	11	5230907	NA	OR51V1 - HBB	283111	3043		        8.98	       15.79	rs7120391-C	rs7120391	0	7120391	Intergenic	1	0.12	5E-9	8.301029995663981	(MCHC)	    .30	[0.20-0.40] unit increase	Illumina [907,954]	N


#FIELDs
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

#VALUE
#0 11/13/2013
#1 23696099
#2 Ding K
#3 05/20/2013
#4 G3 (Bethesda)
#5 http://www.ncbi.nlm.nih.gov/pubmed/23696099
#6 Genetic variants that confer resistance to malaria are associated with red blood cell traits in African-Americans: an electronic medical record-based genome-wide association study.
#7 Red blood cell traits
#8 1904 Afican American individuals
#9 411 Afican American individuals
#10 11p15.4
#11 11
#12 5230907
#13 NA
#14 OR51V1 - HBB
#15 283111
#16 3043	
#17
#18 8.98
#19 15.79
#20 rs7120391-C
#21 rs7120391
#22 0
#23 7120391
#24 Intergenic
#25 1
#26 0.12
#27 5E-9
#28 8.301029995663981
#29 (MCHC)
#30 .30
#31 [0.20-0.40] unit increase
#32 Illumina [907,954]
#33 N

#ATM I'll pick the following field numbers:
#11	12 21 (-)REF     (-)ALT     (-)QUAL    FILTER(7)  INFO(6) PVAL(27)

my $PROBABILITY = 0.00000005; 
my $MINSNPS = 5;

my $infile;
GetOptions(
        'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_GWAScatalog_to_vcf.pl -i=<INFILE>\n";
     print "<INFILE> GWAS catalog TSV file\n";
     print "NOTE: I will only consider SNPs having p < $PROBABILITY and traits with at least $MINSNPS SNPs\n";
     print "CHANGE if needed\n";
     exit 1;
}


#print header
print "##fileformat=VCFv4.1\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");

#there will be duplicates (same all, different pvals)
#need a hash
my %textline;
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^Date/); #header
	#vcf fields
	my $chrom; #11
	my $pos;   #12
	my $id;    #21
	my $missing = '-'; #REF     ALT     QUAL    
	my $filter; #7
	my $info;   #6
	my $pval;
	
	($info, $filter, $chrom, $pos, $id, $pval) = (split /\t/)[6,7,11,12,21,27];
	next if(!$chrom);
	next if($chrom eq '');
	
	#ONLY KEEP GENOME WIDE SIGNIFICANT SNPs
	next unless($pval);
	next if($pval > $PROBABILITY);
	$info = $info . '(p=' . $pval . ')';
	
	#chromosome is in 1...23 format, adjust
	#TODO CURRENTLY PRODUCES b37 style results
	
	if($chrom =~ /23/){
		$chrom = 'X';
	}elsif($chrom =~ /24/){
		$chrom = 'Y';
	}elsif($chrom =~ /25/){
		$chrom = 'M';
	}elsif($chrom =~ /\d+/){
		$chrom = $chrom;
	}else{
		print "Chromosome string $chrom is not recognised. Skipping entry..\n";
		next;
	}
	my $line =  $chrom . "\t" . $pos . "\t" . $id . "\t" . $missing . "\t" . $missing .  "\t" . $missing . "\t" . $filter . "\t" . $info . "\n";
	$textline{$line} = 1; 
}
close $instream;

#2nd pass
#build phenotype->snp-number map based on phenotype
my %phenotype_map;
foreach my $entry (keys %textline){
	my @fields = split (/\t/, $entry);
	if(!$phenotype_map{lc($fields[6])}){
		$phenotype_map{lc($fields[6])} = 1;
	}else{
		$phenotype_map{lc($fields[6])} += 1;		
	}
}

foreach my $item (keys %textline){
	my @fields = split (/\t/, $item);
	next if($phenotype_map{lc($fields[6])} < $MINSNPS); 
	print $item; 
}