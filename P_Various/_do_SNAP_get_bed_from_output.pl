#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#I'm using this Broad program to retrieve snps in LD with a list of snp I input.
#I'm using for example GWAS snps
#http://www.broadinstitute.org/mpg/snap/ldsearch.php
#This outputs a list like the following 
#SNP     Proxy   Distance        RSquared        DPrime  Arrays  Chromosome      Coordinate_HG18
#rs2038256       rs2038256       0       1.000   1.000   A6,I3,I5,I6,I6Q,IM,IMD,IC,ICQ,CYT,IWQ,OE,E1,E11,O54,O5E,OEE,AAE,AAH     chr14   28319100
#rs2038256       rs1950188       1067    1.000   1.000   None    chr14   28320167

#where 
#0 SNP is my query disease-associated snp
#1 Proxy is the SNP in LD
#2 distance is the distance in bp between the two?
#6 chrom
#7 coord_hg18 (remember to liftover)

#bed format is like this
#chrom start end name score

#output will be like this:
#chrom start end rsid_of_snp distance_of_this_snp_from_query

my $infile;
my $infile_gwas; #full gwas vcf - I need this to get coordinates/chr for those snps unrecognised by snap
GetOptions(
        'i1=s'        =>\$infile,
        'i2=s'        =>\$infile_gwas
);
if(!$infile){
	print "USAGE: do_SNAP_get_bed_from_output.pl -i1=<INFILE1> -i2=<INFILE2>\n";
    print "<INFILE1> file .txt from SNAP\n";
    print "<INFILE2> reference GWAS .vcf file\n";
    exit 1;
}
if(!$infile_gwas){
	print "USAGE: do_SNAP_get_bed_from_output.pl -i1=<INFILE1> -i2=<INFILE2>\n";
    print "<INFILE1> file .txt from SNAP\n";
    print "<INFILE2> reference GWAS .vcf file\n";
    exit 1;
}
my $PATH_LIFTOVER = '/net/isi-backup/giuseppe/scripts/liftOver';
#liftOver oldFile map.chain newFile unMapped
my $PATH_CHAIN = '/net/isi-scratch/giuseppe/indexes/Hsap/hg18ToHg19.over.chain';
my ($filename, $directory) = fileparse($infile);
$filename =~ s/(.*)\..*/$1/;
my $outfile_hg18 = $directory . "\/" .  $filename . '_hg18.bed';
my $outfile_hg19 = $directory . "\/" .  $filename . '_hg19.bed';
my $outfile_unmapped = $directory . "\/" .  $filename . '_unmapped.bed';


#build a hash from the vcf catalog with rsID, chromosome, position
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
my %gwas_dictionary;
open (my $instream,  q{<}, $infile_gwas) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#+/);
	my ($chrom, $pos, $id) = (split /\t/)[0,1,2];
	#remove any white spaces from both end of RSid
	$id =~ s/^\s+|\s+$//g;
	$gwas_dictionary{$id} = $chrom . '-' . $pos;
}
close $instream;

open ($instream,  q{<}, $infile) or die("Unable to open $infile : $!");
open (my $hg18_stream,  q{>}, $outfile_hg18) or die("Unable to open $outfile_hg18 : $!");
while(<$instream>){
	$_ =~ s/\R//g; #remove line feed and carriage return
	next if($_ eq '');
	next if($_ =~ /^SNP/); #header
	#if($_ =~ /WARNING/){ print $_, "\n"; next; } #WARNING Query snp not in 1000GenomesPilot1  or WARNING No matching proxy snps found
	
	my ($snp, $snp_proxy, $distance, $rsquared, $chr, $coord) = (split /\t/)[0,1,2,3,6,7];
	if($snp_proxy =~ /^WARNING/){#WARNING Query snp not in 1000GenomesPilot1  or WARNING No matching proxy snps found
		print $_, "\n";
		print "No snps in LD found, storing original snp.\n";
		if(!$gwas_dictionary{$snp}){
			print $snp;
			exit;
		}
		($chr, $coord) = split("-", $gwas_dictionary{$snp});
		print $hg18_stream $chr . "\t" . ($coord-1) . "\t" . $coord . "\t" . $snp . "\t" . '0' . "\n";
		next;
	}
	next if($chr =~ /N\/A/);
	next if(!$coord =~ /N\/A/);
	my $name = $snp . '_' . $snp_proxy; 
	print $hg18_stream $chr . "\t" . ($coord-1) . "\t" . $coord . "\t" . $name . "\t" . $distance . "\n";
}
close $instream;
close $hg18_stream;
#sort file and liftover
system "sort -k1,1V -k2,2g $outfile_hg18 | $PATH_LIFTOVER stdin $PATH_CHAIN $outfile_hg19 $outfile_unmapped";
unlink $outfile_unmapped unless(-s $outfile_unmapped);
