#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#27-1-2014
#This is used to pick all the .err and out  files from the NA<ID> directories produced by do_vcf2diploid_cluster.sh and move them all in a single d_logs directory

my $dir;
my $resultdir = 'd_logs';

GetOptions(
	'd=s'  =>\$dir,
);
if(!$dir){
     print "USAGE:  do_vcf2diploid_group_logs.pl -d=<DIRECTORY>\n";
     print "<DIRECTORY> where the NA* directories are\n";
     exit 1;
}
chdir $dir;
system "mkdir $resultdir";
my $abs_resultdir = $dir . '/' . $resultdir;
my @dirs = <NA*>;

foreach my $in_dir (@dirs) {
	print $in_dir, "\n";
	#system "cd $in_dir";
	my $abs_in_dir = $dir . '/' . $in_dir;
	chdir $abs_in_dir;
	system "cp *.err $abs_resultdir";
	system "cp *.out $abs_resultdir";
}
