use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#INFO 30/1/2015
#This converts the GWAS catalog .txt table found here
#http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz
#in a bed to be used with bedtools (I want to intersect it with my VDR-BVs to build forest plots for FIGURE 4 of the paper)
#in a minimal vcf to be used with bedtools 
#The MINIMUM vcf understood by bedtools contains the following fields
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#So create these

#SAMPLE INPUT
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

#ATM I'll pick the following field numbers:
#1 3 4 (-)REF     (-)ALT     (-)QUAL  FILTER(10)  INFO(9) PVAL(17)

my $PROBABILITY = 0.00000005; 
#my $MINSNPS = 5;

my $infile;
GetOptions(
        'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_GWAScatalog_to_vcf.pl -i=<INFILE>\n";
     print "<INFILE> GWAS catalog UCSC.gz file\n";
     #print "NOTE: I will only consider SNPs having p < $PROBABILITY and traits with at least $MINSNPS SNPs\n";
     print "NOTE: I will only consider SNPs having p < $PROBABILITY\n";
     print "CHANGE if needed\n";
     exit 1;
}


#print header
print "##fileformat=VCFv4.1\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
#there will be duplicates (same all, different pvals)
#need a hash
my %textline;

tie *FILE,   'IO::Zlib', $infile, "rb";
while (<FILE>)	{ 
#open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
#while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^Date/); #header
	#vcf fields
	my $chrom; #1
	my $pos;   #3
	my $id;    #4
	my $missing = '-'; #REF     ALT     QUAL    
	my $filter; #10
	my $info;   #9
	my $pval; #17
	
	#ATM I'll pick the following field numbers:
#1 3 4 (-)REF     (-)ALT     (-)QUAL  FILTER(10)  INFO(9) PVAL(17)
	($chrom, $pos, $id, $info, $filter, $pval) = (split /\t/)[1,3,4,9,10,17];
	next if(!$chrom);
	next if($chrom eq '');
	
	$chrom =~ s/chr(.*)/$1/; #ucsc is hg19, I want b37 output
	
	#ONLY KEEP GENOME WIDE SIGNIFICANT SNPs
	next unless($pval);
	next if($pval =~ /NS/);
	next if($pval eq 'E');
	
	next if($pval > $PROBABILITY);
	$info = $info . '(p=' . $pval . ')';
	
	my $line =  $chrom . "\t" . $pos . "\t" . $id . "\t" . $missing . "\t" . $missing .  "\t" . $missing . "\t" . $filter . "\t" . $info . "\n";
	$textline{$line} = 1; 
}
#close $instream;
close FILE;

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
	#next if($phenotype_map{lc($fields[6])} < $MINSNPS); 
	print $item; 
}