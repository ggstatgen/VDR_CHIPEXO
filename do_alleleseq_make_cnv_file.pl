#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#create the file shown here
#http://info.gersteinlab.org/AlleleSeq
#(d) CNV file; a set of normalized read depth values for all snp locations from a separate genomic sequencing experiment. This reports the read depth at that snp #compared to overall coverage. This is used to filter out locations with very low or high coverage, which would tend to indicate copy number variation. The file should #be in this format (with header):
#
#chrm    snppos  rd
#1       52066   0.902113
#1       695745  0.909802
#1       742429  0.976435
#
#This can be generated from CNVnator (Abyzov et al. 2011). 

my $infile_cnv;
my $infile_snp;
my %outdata;
my $BEDTOOLS  = '/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin';

GetOptions(
	"cnv=s" => \$infile_cnv,
	'vcf=s'  =>\$infile_snp,
);
if(!$infile_cnv){
     print "USAGE: do_alleleseq_make_cnv_file.pl -cnv=<INFILE1> -vcf=<INFILE2>\n";
     print "<INFILE1> CNVnator tsv output\n";
     print "<INFILE2> vcf/vcf.gz file with snps for the trio that includes the child of interest\n";
     exit 1;
}
if(!$infile_snp){
     print "USAGE: do_alleleseq_make_cnv_file.pl -cnv=<INFILE1> -vcf=<INFILE2>\n";
     print "<INFILE1> CNVnator tsv output\n";
     print "<INFILE2> vcf/vcf.gz file with snps for the trio that includes the child of interest\n";
     exit 1;
}

#$infile_snp="/net/isi-scratch/giuseppe/VDR/VARIANTS/omni25_hg19_b37/d_TRIO_DATA/Omni25_genotypes_2141_samples.hg19_OMNI_25_GEN_PASS_YRI_Y111.vcf.gz";
#format: 
#CNV_type        coordinates     CNV_size        normalized_RD   p-val1  p-val2  p-val3  p-val4  q0
#deletion        chr1:1-10000    10000   0       1.59373e-11     1.08391e-55     1.99216e-11     8.31e-43        -1
#$infile_cnv="/net/isi-scratch/giuseppe/VDR/1000g_GENOME_ALIGNMENTS/_NA19189/cnvnator_NA19189.mapped_bin500_calls.tsv";

#my ($filename, $directory) = fileparse($infile);
#$filename =~ s/(.*)\..*/$1/;
#my $outfile = $directory . "\/" .  $filename . '_segments.bed';

#------------
#1 - create bed from cnv file
#------------
open (my $instream,     q{<}, $infile_cnv) or die("Unable to open $infile_cnv : $!");
my ($filename, $directory) = fileparse($infile_cnv);
$filename =~ s/(.*)\..*/$1/;
my $outfile_cnv_bed          = $directory . "\/"                   . $filename . '.bed';
my $outfile_intersect_bed    = $directory . "\/" . 'intersect_'    . $filename . '.bed';
my $outfile_no_intersect_bed = $directory . "\/" . 'no_intersect_' . $filename . '.bed';

open (my $temp_stream,  q{>}, $outfile_cnv_bed) or die("Unable to open $outfile_cnv_bed : $!");
my $counter_del = 1;
my $counter_dup = 1;
while(<$instream>){
	next if($_ =~ /^CNV_type/);	#header
	my $chr; my $start; my $end; my $name;
	chomp;
	my ($type, $coords, $norm_rd) = (split /\t/)[0,1,3];
	#create bed chrom start end name score
	if ($coords =~ /(.+):(\d+)-(\d+)/){
		$chr = $1; $start = $2; $end = $3;
	}else{
		print "$coords: format not recognised. Skipping.\n";
		next;
	}
	if($type =~ /deletion/){
		$name = $type . '_' . $counter_del;
		$counter_del++;
	}elsif($type =~ /duplication/){
		$name = $type . '_' . $counter_dup;
		$counter_dup++;
	}else{
		print "Unknown type: $type. Skipping..\n";
		next;		
	}
	$chr =~ s/chr(.+)/$1/;
	my $line = $chr . "\t" . $start . "\t" . $end . "\t" . $name . "\t" . $norm_rd; 
	print $temp_stream $line, "\n";
}
close $instream;
close $temp_stream;
#------------
#2 - intersect
#------------
#bedtools -wo -a -b
system "sort -k1,1V -k2,2g $outfile_cnv_bed | $BEDTOOLS/bedtools intersect -wo -a $infile_snp -b stdin > $outfile_intersect_bed";
#bedtools -v -a -b
system "sort -k1,1V -k2,2g $outfile_cnv_bed | $BEDTOOLS/bedtools intersect -v  -a $infile_snp -b stdin > $outfile_no_intersect_bed";
#-----------
#3 - go through the two files and fill hash
#-----------
#output must be in this format:
#chrm    snppos  rd
#1       52066   0.902113
#1       695745  0.909802
#1       742429  0.976435
open ($instream,     q{<}, $outfile_intersect_bed) or die("Unable to open $outfile_intersect_bed : $!");
while(<$instream>){
	chomp;
	my ($chr, $snppos, $rd ) = (split /\t/)[0,1,14];
	$chr =~ s/chr(.+)/$1/;
	my $line = $chr . "\t" . $snppos;
	$outdata{$line} = $rd;
}
close $instream;

open ($instream,     q{<}, $outfile_no_intersect_bed) or die("Unable to open $outfile_no_intersect_bed : $!");
while(<$instream>){
	chomp;
	my ($chr, $snppos) = (split /\t/)[0,1];
	$chr =~ s/chr(.+)/$1/;
	my $line = $chr . "\t" . $snppos;
	$outdata{$line} = 1;
}
close $instream;

print "chrm\tsnppos\trd\n";
foreach my $item (sort keys %outdata){
	print $item . "\t" . $outdata{$item} . "\n";
}
unlink $outfile_no_intersect_bed;
unlink $outfile_intersect_bed;
unlink $outfile_cnv_bed;
