#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use List::Util qw(min max);


#I carry out two peak calls, one with d=12 one with d=24
#I want to merge these peaksets as described in the paper, so that:

#1 if a peak is present only in 1 bed, pick it
#2 if a peak is present in both, intersect it

my $infile1;
my $infile2;
my $temp1;
my $temp2;
my $temp3;
my $temp4;

GetOptions(
        'i1=s'  => \$infile1,
        'i2=s'  => \$infile2
);
if(!$infile1){
     print "USAGE: do_macs2_merge_peaksets.pl -i1=<INFILE1> -i2=<INFILE2>\n";
     print "<INFILE1> bed file 1\n";
     print "<INFILE2> bed file 2\n";
     print "NOTE: INFILE1 interval will be reported when overlap with INFILE2 interval is present!!\n";
     exit 1;
}
if(!$infile2){
     print "USAGE: do_macs2_merge_peaksets.pl -i1=<INFILE1> -i2=<INFILE2>\n";
     print "<INFILE1> bed file 1\n";
     print "<INFILE2> bed file 2\n";
     print "NOTE: INFILE1 interval will be reported when overlap with INFILE2 interval is present!!\n";
     exit 1;
}

my $BEDTOOLS_PATH  = '/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin/';

#%output goes in first dir
my($basename, $directory) = fileparse($infile1);
$basename =~ s/(.*)\..*/$1/;
my $ID;
if($basename =~ /.*(GM\d{5}).*/){
	$ID = $1;
}else{
	$ID = $basename;
}
$temp1 = $directory . $ID  . '_temp1.bed';
$temp2 = $directory . $ID  . '_temp2.bed';
$temp3 = $directory . $ID  . '_temp3.bed';
$temp4 = $directory . $ID  . '_temp4.bed';


#get those peaks in 1 which DO NOT overlap with 2
#get those peaks in 2 which DO NOT overlap with 1
system "$BEDTOOLS_PATH/bedtools intersect -v -a $infile1 -b $infile2 > $temp1";
system "$BEDTOOLS_PATH/bedtools intersect -v -a $infile2 -b $infile1 > $temp2";

#same result??
#system "$BEDTOOLS_PATH/bedtools subtract -A -a $infile1 -b $infile2 > $temp1";
#system "$BEDTOOLS_PATH/bedtools subtract -A -a $infile2 -b $infile1 > $temp1";

#merge and sort
system "cat $temp1 $temp2 | sort -k1,1V -k2,2g > $temp3";

#intersect initial
#system "$BEDTOOLS_PATH/bedtools intersect  -a $infile1 -b $infile2 > $temp4";
system "$BEDTOOLS_PATH/bedtools intersect -wa -u -a $infile1 -b $infile2 > $temp4";


#merge and sort
system "cat $temp3 $temp4 | sort -k1,1V -k2,2g";
#result on stdout

unlink $temp1; 
unlink $temp2; 
unlink $temp3; 
unlink $temp4; 
