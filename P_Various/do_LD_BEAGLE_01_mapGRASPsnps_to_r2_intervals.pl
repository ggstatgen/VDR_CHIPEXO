#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#I got a set of LD intervals based on a r^2 definition from here
#https://data.broadinstitute.org/srlab/BEAGLE/1kG-beagle-release3/ld_intervals/

#see also
#http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5/READ_ME_beagle_ref
#format is
#chr1    10583   10583   rs58108140
#chr1    13957   13957   1:13957

#The author says (mail:slowikow@broadinstitute.org )
# I take each SNP in the 1000G project and find the two most distant SNPs upstream and downstream with r2 greater than some threshold. The second and third columns are the positions of the two most distant SNPs (within 1Mb).

#So I want to map my GRASP genome wide snps to their intervals and rewrite in a bed file
#only snps with p < 5e-8; only traits with 5+ snps at least
my $INPUT_GRASP = "/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/GRASP2final.vcf.gz";
#These are 23,555 unique snps in vcf

#output:
#a bed file similar to the grasp vcf
#chrom ld_start ld_stop  name score
#name: disease
#score: pval of the index snp

#method. Create a  hash:
#key: snp; value: interval
#or
#key
#chr:pos: value: interval
#query both
#if the interval around a snp is 1nt, you don't need to save it into the hash

my $infile;
GetOptions(
        'i=s'      =>\$infile
); 
if(!$infile){
	print "USAGE: do_LD_mapGRASPsnps_to_BROADrsquare_intervals.pl -i=<INFILE>\n";
    print "<INFILE> bed.gz file of intervals from the BROAD\n";
    print "eg /net/isi-mirror/1000_Genomes/POPULATION_PANELS/BEAGLE_r2_BLOCKS/TGP2011.0_8.EUR_hg19.bed.gz\n";
    exit 1;
}

my %LD_dataset;
#chr > rsid/pos > interval

tie *FILE,   'IO::Zlib', $infile, "rb";
while (<FILE>)	{ 
	chomp;
	next if($_ eq '');
	next if($_ =~ /^#/);

	my @fields = split("\t", $_);
	#Interval: if the interval is 1nt long, it means it has nothing else in ld. I may skip it.
	#I will not find  the corresponding SNP in the GRASP2 when I query the hash, I will save its coordinate
	next if($fields[2] eq $fields[1]);
	
	my $interval = $fields[1] . '-' . $fields[2];	
	#clean chromosome
	my $chr = $fields[0];
	$chr =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19
	my $id = $fields[3]; #this can be an rsid or a 1:pos
	$LD_dataset{$chr}{$id} = $interval;
}
close FILE;



#now get the rsids and intervals in 1:pos format, query the hash and get a bed
##fileformat=VCFv4.1
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#1       807512  rs10751454      -       -       -       Cardiovascular disease prevalence       Cardiovascular disease (CVD);Myocardial infarction (MI);Mortality;Cancer(p=2.20E-13)
#1       807512  rs10751454      -       -       -       Systolic blood pressure (SBP)   Cardiovascular disease (CVD);Myocardial infarction (MI);Mortality;Cancer(p=5.40E-21)
#1       807512  rs10751454      -       -       -       Total cholesterol       Cardiovascular disease (CVD);Myocardial infarction (MI);Mortality;Cancer(p=7.70E-09)
tie *FILE,   'IO::Zlib', $INPUT_GRASP, "rb";
while (<FILE>)	{
	chomp;
	next if($_ eq '');
	next if($_ =~ /^#/);
	
	my @fields = split("\t", $_);
	#create alternative id
	my $alt_id = $fields[0] . ':' . $fields[1];
	
	if($LD_dataset{$fields[0]}{$fields[2]}){
		my @ld_boundaries = split('-', $LD_dataset{$fields[0]}{$fields[2]});
		print $fields[0] . "\t" . ($ld_boundaries[0]-1) . "\t" . $ld_boundaries[1] . "\t" . $fields[6] . "\t" . $fields[7] . ';' . $fields[2] . "\n"; 
	}elsif($LD_dataset{$fields[0]}{$alt_id}){
		my @ld_boundaries = split('-', $LD_dataset{$fields[0]}{$alt_id});
		print $fields[0] . "\t" . ($ld_boundaries[0]-1) . "\t" . $ld_boundaries[1] . "\t" . $fields[6] . "\t" . $fields[7] . ';' . $fields[2] . "\n"; 		
	}else{
		print $fields[0] . "\t" . ($fields[1]-1) . "\t" . $fields[1] . "\t" . $fields[6] . "\t" . $fields[7] . ';' . $fields[2] .  "\n";
	}		
}
close FILE;