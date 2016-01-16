#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#this tools wants signal in a sgl file
#http://info.gersteinlab.org/ACT_Tool
#I have wigs out macs (or also bedgraphs or bigwigs)
#
#sgl looks like this:
#chr22   15241916        8
#chr22   15241935        9
#chr22   15241938        11
#chr22   15241939        12
#chr22   15241964        11
#chr22   15241986        10
#chr22   15242006        9
#chr22   15242013        7
#chr22   15242048        8
#chr22   15242050        10
#chr22   15242066        9
#chr22   15242079        10
#chr22   15242082        9
#chr22   15242116        8
#chr22   15242131        9
#chr22   15242135        8
#chr22   15242138        6
#chr22   15242139        5
#
#wig is like this
#track type=wiggle_0 name="/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_CEU/MACS14_VDR_CEU_merged.sorted_treat_all" description="Extended tag pileup from MACS version 1.4.2 20120305 for every 10 bp"
#variableStep chrom=chr1 span=10
#10001   3
#10011   6
#10021   8
#10031   5
#10231   1
#10241   5
#10251   9
#10261   8
#10271   3
#10281   2
#10321   2
#10331   4
#10341   2
#


my $infile;
GetOptions(
        'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_wig_to_sgl.pl -i=<FILE>\n";
     print "<FILE> input wig file\n";
     exit 1;
}

my $chr; my $span;
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	next if($_ =~ /^track/);
	if($_ =~ /^variableStep chrom=(chr.+) span=(\d+)/){
		$chr = $1;
		$span = $2;#not used
		next;
	}
	my ($coord, $score) = (split /\t/)[0, 1];
	print "$chr\t$coord\t$score";
}
close $instream;

#my $chr="chr";
#my $start=-1;
#my $sig=0;
#while(<$instream>)
#{
#	my @line = split/\s+/, $_;
#	my $temp_chr = $line[0];
#	my $temp_start = $line[1]-1;
#	my $temp_sig = $line[2];	
#	if ($chr ne $temp_chr){
#		$chr = $temp_chr;
#		$start = $temp_start + 1;
#		$sig = $temp_sig;
#	}
#	else{
#		print "$chr\t$start\t$temp_start\t$sig\n";
#		$start = $temp_start + 1;
#		$sig = $temp_sig;
#	}
#	
#}
#close $instream;
