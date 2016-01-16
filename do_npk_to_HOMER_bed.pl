#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#10/10/2013
#I want to use the genome annotation tools in HOMER
##http://biowhat.ucsd.edu/homer/ngs/annotation.html
##the guide says the program requires a bed formatted as follows:
#BED files should have at minimum 6 columns (separated by TABs, additional columns will be ignored)
#
#Column1: chromosome
#Column2: starting position
#Column3: ending position
#Column4: Unique Peak ID
#Column5: not used
#Column6: Strand (+/- or 0/1, where 0="+", 1="-")
#
#In theory, HOMER will accept BED files with only 4 columns (+/- in the 4th column), and files without unique IDs, but this is NOT recommended.
#
#As MACS doesn't give strand info in the peaks, I will be producing temporary beds with four fields and see how it goes
#
#The author was emailed and said
#
#I would either fake the strand (just add + in the 4th or 6th column), or just give it a 3 column bed file (chr start end).  It should work with either - for chip-seq the strand isn't really important, but it's good to have to maintain compatibility with things and for HOMER sometimes you may want to use a peak/BED that does have strand information, such as a transcript.
#
#
##This file turns a bed produced with MACS2 as follows
#chr1    974244  974310  /net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_CEU/VDR_CEU_merged.sorted_peak_1   173
#chr1    1000121 1000179 /net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_CEU/VDR_CEU_merged.sorted_peak_2   124
#chr1    1080874 1080936 /net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_CEU/VDR_CEU_merged.sorted_peak_3   180
#chr1    1209207 1209369 /net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_CEU/VDR_CEU_merged.sorted_peak_4   443
#
#into
#chr1    974244  974310  /net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_CEU/VDR_CEU_merged.sorted_peak_1	.	+
#
#



my $infile;
GetOptions(
	'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_npk_to_bed.pl -i=<FILE>\n";
     print "<FILE> input narrowPeak file\n";
     exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");

while(<$instream>){
	my $peak_ID;
	my ($chrom, $chromStart, $chromEnd,$peak_name) = (split /\t/)[0, 1, 2,3];
	if($peak_name =~ /.*(peak_\d+)$/){
		$peak_ID = $1;
	}else{
		$peak_ID = $peak_name;
	}
	my $newline = $chrom . "\t" . $chromStart . "\t" . $chromEnd . "\t" . $peak_ID . "\t" . "" . "\t" . "+" ; 
	print $newline, "\n";
}
close $instream;
