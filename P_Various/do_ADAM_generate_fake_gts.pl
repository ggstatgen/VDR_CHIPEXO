#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Zlib;

#28/9/2015
#Adam wants me to create a reference where all alternate alleles are mapped in. I can use vcf2diploid
##fileformat=VCFv4.0
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele, ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/technical/reference/ancestral_alignments/README">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=HM2,Number=0,Type=Flag,Description="HapMap2 membership">
##INFO=<ID=HM3,Number=0,Type=Flag,Description="HapMap3 membership">
##reference=human_b36_both.fasta
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12891 NA12892 NA12878
#20      9098    rs6078030       C       T       .       PASS    AA=.;DP=108     GT:GQ:DP        0|0:81:27       1|1:73:29       0|1:100:48

#for it to work, I need genotypes, which should have the paternal and maternal.
#I think in order to create a full alternate, I could simply generate a fake genotype column with a 1|0 or a 0|0

#java -jar ../vcf2diploid.jar \
#          -id NA12878 \
#          -chr human_chr20_hg18.fa \
#	       -vcf CEU.trio.chr20.2010_03.genotypes.vcf

#INPUT: vcf
#OUTPUT: vcf with fake gt line

my $infile;
my $SNPSONLY;
GetOptions(
        'i=s'      =>\$infile,
        'snpsonly' =>\$SNPSONLY
);

#$infile = '/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/ALLELESEQ_20samples_merged.vcf.gz'; #b37
#$SNPSONLY = '1';

if(!$infile){
     print "USAGE: do_ADAM_generate_fake_gts.pl -i=<INFILE> -snpsonly\n";
     print "<INFILE> 1kg .vcf.gz file\n";
     print "(optional) <snpsonly> do not include cnvs\n";
     exit 1;
}

tie *FILE,   'IO::Zlib', $infile, "rb";
while (<FILE>)	{ 
	if($_ =~ /^\#\#/){
		print $_;
		next;
	}
	
	chomp;
	if($_ =~ /^\#CHROM/){
		#add a fake sample
		print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_c\tSAMPLE_p\tSAMPLE_m\n";
		next;
	}
	#my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = (split /\t/)[0,1,2,3,4,5,6,7,8];
	my @data = (split /\t/)[0,1,2,3,4,5,6,7,8];
	
	if(!$data[8]){
		$data[8] = 'GT';
	}
	
	if($SNPSONLY){
		next if(length($data[3]) > 1);
		next if(length($data[4]) > 1);
	}

	#some filtering, change based on what Adam thinks
	#next unless ($data[2] =~ /^rs/i);
	next unless ($data[6] =~ /PASS/);
	
	print STDOUT join("\t", @data) . "\t" . '1|0' . "\t" . '1|1' . "\t" . '0|0' . "\n";
}
close FILE;