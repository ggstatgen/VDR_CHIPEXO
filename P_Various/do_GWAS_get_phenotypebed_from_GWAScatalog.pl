#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $TOOL_PATH       = '/net/isi-scratch/giuseppe/tools';
my $BEDTOOLS  = $TOOL_PATH . '/bedtools-2.17.0/bin';
my $CHROM_SIZE_PATH = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes';

#INFO
#I want to extract a bed file with all the intervals around all the snps for a particular phenotype (for example "Multiple Sclerosis")
#the interval around the snps must be selectable
#I need this "one disease" data to intersect it with my peaks
#INPUT: gwas catalog (ideally, one with added MS data)
#output 1: bed with all the disease snps
#output 2: bed with all the disease snps and window around them

#this should be very similar to do_GWAScatalog_to_bed.pl

my $infile;
my $disease;
my $SLOP;
GetOptions(
        'i=s'      =>\$infile,
        'd=s'      =>\$disease,
        'slop=i'   =>\$SLOP
);
if(!$infile){
     print "USAGE: do_GWAS_get_phenotypebed_from_GWAScatalog.pl -i=<INFILE> -d=<DISEASE> -slop=<BASEPAIRS>\n";
     print "<INFILE> GWAScatalog txt file\n";
     print "<DISEASE> Disease/Trait string (e.g. multiple sclerosis)\n";
     print "<BASEPAIRS> number of bp requested upstream (downstream) of each SNP\n";
     exit 1;
}
if(!$disease){
     print "USAGE: do_GWAS_get_phenotypebed_from_GWAScatalog.pl -i=<INFILE> -d=<DISEASE> -slop=<BASEPAIRS>\n";
     print "<INFILE> GWAScatalog txt file\n";
     print "<DISEASE> Disease/Trait string (e.g. multiple sclerosis)\n";
     print "<BASEPAIRS> number of bp requested upstream (downstream) of each SNP\n";
     exit 1;
}
if(!$SLOP){
     print "USAGE: do_GWAS_get_phenotypebed_from_GWAScatalog.pl -i=<INFILE> -d=<DISEASE> -slop=<BASEPAIRS>\n";
     print "<INFILE> GWAScatalog txt file\n";
     print "<DISEASE> Disease/Trait string (e.g. multiple sclerosis)\n";
     print "<BASEPAIRS> number of bp requested upstream (downstream) of each SNP\n";
     exit 1;
}
my ($filename, $directory) = fileparse($infile);
$filename =~ s/(.*)\..*/$1/;
my $disease_string = $disease;
$disease_string =~ s/ /_/g;
my $outfile = $directory . "\/" .  $filename . '_' . $disease_string . '.bed';
my $outfile_slop = $directory . "\/" .  $filename . '_' . $disease_string . '_' . $SLOP . 'kb.bed';

#open the input and select only the entries with the required disease
#collect unique entries (study,disease/trait,chromosome,position,snp id) in hash
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
my %unique_entries;
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^Date/); #header
	
	#There are duplicate SNP-phenotype lines.
	my ($dis_trait, $chrom, $pos, $snp_id) = (split /\t/)[7,11,12,21];
	#maybe we need some more refined regex here
	next unless(lc($dis_trait) =~ /\Q$disease/);
	next if($snp_id eq '');
	next if(!$chrom);
	next if($chrom eq '');
	
	my $line = $chrom . "\t" . $pos . "\t" . $snp_id . "\t" . $dis_trait;
	if(!$unique_entries{$line}){
		$unique_entries{$line} = 1;
	}else{
		$unique_entries{$line} += 1;		
	}
}
close $instream;

#now you can create the first bed, with only the snp positions
open (my $out_stream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
foreach my $item (sort keys %unique_entries){
	my ($chrom, $pos, $snp_id, $dis_trait) = split (/\t/, $item);
	if($chrom =~ /23/){
		$chrom = 'chrX';
	}elsif($chrom =~ /24/){
		$chrom = 'chrY';
	}elsif($chrom =~ /25/){
		$chrom = 'chrM';
	}elsif($chrom =~ /\d+/){
		$chrom = 'chr' . $chrom;
	}else{
		print "Chromosome string $chrom is not recognised. Skipping entry..\n";
		next;
	}
	#create bed name ( snp_id (disease/trait)   )
	my $bed_name = $snp_id . '(' . $dis_trait . ')';
	#use 0-based start coordinates, which are the proper BED format (pos-1, pos)
	my $line = $chrom . "\t" . ($pos-1)  . "\t" . $pos . "\t" . $bed_name . "\n";
	print $out_stream $line;
}
close $out_stream;

#now extend each snp by $SLOP using bedtools, merge overlapping intervals, save in second file
system "$BEDTOOLS/bedtools slop -b $SLOP -i $outfile -g $CHROM_SIZE_PATH |  sort -k1,1V -k2,2g | $BEDTOOLS/bedtools merge -n -i stdin > $outfile_slop";
#open ($out_stream,  q{>}, $outfile_slop) or die("Unable to open $outfile_slop : $!");
#close $out_stream;
