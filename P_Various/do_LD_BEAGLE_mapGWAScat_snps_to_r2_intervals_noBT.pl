#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#like the other but this does not use bedtools merge so it should be faster. However I trust bedtools merge more.

#Figure 6 - gwas effect. 
#I want to show results not only for the VDR-hcBVs that are GWAScat SNPs, but also for those which are in perfect LD (r^2 = 1) with these
#This script will map GWAScatalog SNPs to their perfect LD blocks

#needs to be improved if you want to pick traits with at least 5 snps

#I got a set of LD intervals based on a r^2 definition from here
#https://data.broadinstitute.org/srlab/BEAGLE/1kG-beagle-release3/ld_intervals/

#see also
#http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5/READ_ME_beagle_ref
#format is
#chr1    10583   10583   rs58108140
#chr1    13957   13957   1:13957

#The author says (mail:slowikow@broadinstitute.org )
# I take each SNP in the 1000G project and find the two most distant SNPs upstream and downstream with r2 greater than some threshold. The second and third columns are the positions of the two most distant SNPs (within 1Mb).

#So I want to map my GWAScat genome wide snps to their intervals and rewrite in a bed file
#my $INPUT_GWAScat = "/net/isi-scratch/giuseppe/indexes/GWAS_Gwascat/UCSC_gwas_catalog/gwasCatalog.hg19.UCSC.gz";
#new, includes the Beecham MS variants
my $INPUT_GWAScat = "/net/isi-scratch/giuseppe/indexes/GWAS_Gwascat/UCSC_gwas_catalog/gwasCatalog.hg19.UCSC_plus_Beecham92-MS.gz";
#(filter out non genome wide significant from these)

#output:
#a bed file similar to the output of do_GWAScatalog_UCSC_to_bed.pl
#it NEEDS to contain EVERYTHING (OR,CI, etc)
#YOU NEED TO SAVE the SNP coordinate (hg19) in the 4 column somewhere
#BED FIELDS
#chrom	start	end	name	score
#1	ldstart	ldstop	snpcord_4-22(__ separated)


#method. Create a  hash:
#key: snp; value: interval
#or
#key
#chr:pos: value: interval
#query both
#if the interval around a snp is 1nt, you don't need to save it into the hash

my $PROBABILITY = 0.00000005;

my $infile;
my $outdir;
GetOptions(
        'i=s'      =>\$infile
); 

#$infile = "/net/isi-mirror/1000_Genomes/POPULATION_PANELS/BEAGLE_r2_BLOCKS/TGP2011.1_0.EUR_hg19_clean.chr16.bed.gz";
if(!$infile){
	print "USAGE: do_LD_BEAGLE_mapGWAScat_snps_to_r2_intervals.pl -i=<INFILE> -out=<OUTDIR>\n";
    print "<INFILE> bed.gz file of intervals from the BROAD\n";
    print "<OUT> output directory\n";
    print "eg <INFILE> /net/isi-mirror/1000_Genomes/POPULATION_PANELS/BEAGLE_r2_BLOCKS/TGP2011.0_8.EUR_hg19.bed.gz\n";
    print "GWAS snps will be from $INPUT_GWAScat\n";
    exit 1;
}
my %LD_dataset;
#chr > rsid/pos > interval

#ATTENTION - the same SNP (field 3) MIGHT be reported multiple times with different LD blocks around
#eg
#chr16   85983326        86021505        rs17445836
#chr16   86017663        86017663        rs17445836
#query a hash of SNP presence. If snp already seen, merge, compute union of intervals

tie *FILE,   'IO::Zlib', $infile, "rb";
while (<FILE>)	{ 
	chomp;
	next if($_ eq '');
	next if($_ =~ /^#/);

	my @fields = split("\t", $_);
	#Interval: if the interval is 1nt long, it means it has nothing else in ld. I may skip it.
	#I will not find  the corresponding SNP in the GWAScat when I query the hash, I will save its coordinate
	next if($fields[2] eq $fields[1]);
	my $id = $fields[3]; #this can be an rsid or a 1:pos
	my $chr = $fields[0];
	my $start = $fields[1];
	my $stop = $fields[2];
	my $interval; #which interval to save?	#merging with bedtools atm
	
	unless($LD_dataset{$chr}{$id}){
		$interval = $start . '-' . $stop;
		$LD_dataset{$chr}{$id} = $interval;
		next;
	}
	#if here, the rsID had been stored already with some LD block info. Retrieve and compare with new ld block info
	my $stored_interval = $LD_dataset{$chr}{$id};
	my ($stored_start, $stored_stop) = split('-',$stored_interval);
	
	#disjoint
	if( (($start <  $stored_start) && ($stop < $stored_stop))  ||    
	    (($start >  $stored_start) && ($stop > $stored_stop))  ){
		print STDERR "Stored and new intervals are disjoint, it should not be here. Aborting.\n";
		exit -1;		    
	}
	#we now know that intervals overlap
	my $final_start;
	my $final_stop;
	if($start <= $stored_start){
		$final_start = $start;
		if($stop >= $stored_stop){
			$final_stop = $stop
		}else{
			$final_stop = $stored_stop;
		}
	}else{ #start > stored_start
		$final_start = $stored_start;
		if($stop >= $stored_stop){
			$final_stop = $stop
		}else{
			$final_stop = $stored_stop;
		}
	}
	$interval = $final_start . '-' . $final_stop;
	$LD_dataset{$chr}{$id} = $interval;
}
close FILE;


#now get the rsids and intervals in 1:pos format, query the hash and get a bed
#save all fields from the GWAS catalog
#filter on GWAS genomewide significance
#592     chr1    1005805 1005806 rs3934834       19851299        Johansson A     2009-10-22      Obesity (Silver Spring) Linkage and genome-wide association analysis of obesity-related phenotypes: association of weight with the MGAT1 gene.  Body mass index 1,079 South Tyrolian individuals, 790 Dutch founder individuals, 2,060 European ancestry individuals    NA      1p36.33 NR      rs3934834-G     0.80    6E-7    (females + males)       .11     [NR] kg increase        Illumina [318,237]      N

#snp pos is field 3,save it in the bed

tie *FILE,   'IO::Zlib', $INPUT_GWAScat, "rb";
while (<FILE>)	{
	chomp;
	next if($_ eq '');
	my @data = split("\t", $_);
	
	#check significance of snp; if not genome wide significant we don't want the interval
	next unless($data[17] =~ /\d+/);
	next if($data[17] > $PROBABILITY);
	
	my $unknown = shift @data;
	my $chr = shift @data;
	#$chr =~ s/chr(.*)/$1/;#convert hg19 to b37
	
	my $start = shift @data; 
	my $snp_pos = shift @data;
	
	my $rs_id = $data[0];
	my $alt_id = $chr . ':' . $snp_pos; #remember that the input file is a 0-based bed
	my $gwas_data = join('__', @data);
	$gwas_data = $snp_pos . '__' . $gwas_data; #add the actual snp position at the beginning of the vector of gwas info
	
	if($LD_dataset{$chr}{$rs_id}){
		my @ld_boundaries = split('-', $LD_dataset{$chr}{$rs_id});
		print $chr . "\t" . ($ld_boundaries[0]-1) . "\t" . $ld_boundaries[1] . "\t" .  $gwas_data  . "\n"; 
	}elsif($LD_dataset{$chr}{$alt_id}){
		my @ld_boundaries = split('-', $LD_dataset{$chr}{$alt_id});
		print $chr . "\t" . ($ld_boundaries[0]-1) . "\t" . $ld_boundaries[1] . "\t" .  $gwas_data . "\n"; 		
	}else{
		print $chr . "\t" . $start  . "\t" . $snp_pos . "\t" . $gwas_data .  "\n";
	}	
}
close FILE;