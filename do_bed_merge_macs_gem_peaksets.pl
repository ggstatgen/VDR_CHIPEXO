#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use List::Util qw(min max);


#I carry out two peak calls, one with gem and one with macs (consensus of other peak calls)
#I want to merge these peaksets as described in the paper, so that:

#1 if a peak is present only in 1 bed, pick it
#2 if a peak is present in both, get the union of them
#3 merge close peaks??

my $infile1;
my $infile2;
my $distance;
GetOptions(
        'i1=s'  => \$infile1,
        'i2=s'  => \$infile2,
	'd=i'   => \$distance
);
if(!$infile1){
     print "USAGE: do_bed_merge_macs_gem_peaksets.pl -i1=<INFILE1> -i2=<INFILE2> -d=<DISTANCE>\n";
     print "<INFILE1> bed file 1\n";
     print "<INFILE2> bed file 2\n";
     print "(optional)<DISTANCE> whether to merge features distant <DISTANCE>bp apart. Default: do not merge non-overlapping features\n";
     exit 1;
}
if(!$infile2){
     print "USAGE: do_bed_merge_macs_gem_peaksets.pl -i1=<INFILE1> -i2=<INFILE2> -d=<DISTANCE>\n";
     print "<INFILE1> bed file 1\n";
     print "<INFILE2> bed file 2\n";
     print "(optional)<DISTANCE> whether to merge features distant <DISTANCE>bp apart. Default: do not merge non-overlapping features\n";
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

#I think a cat of the two files plus a merge would work here
if($distance){
	system "cat $infile1 $infile2 | sort -k1,1 -k2,2n | $BEDTOOLS_PATH/bedtools merge -d $distance -i stdin";
}else{
	system "cat $infile1 $infile2 | sort -k1,1 -k2,2n | $BEDTOOLS_PATH/bedtools merge -i stdin";
}
