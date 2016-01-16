#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#17-4-2013
#This is used to pick all the *_1_GEM_events.bed from the _outputs directory produced by GPS and move them all in a single d_beds directory


my $dir;
my $resultdir = 'd_beds';

GetOptions(
	'd=s'  =>\$dir,
);
if(!$dir){
     print "USAGE:  do_GEM_group_beds.pl -d=<DIRECTORY>\n";
     print "<DIRECTORY> where the *_output directories are\n";
     exit 1;
}
chdir $dir;
system "mkdir $resultdir";
my $abs_resultdir = $dir . '/' . $resultdir;
my @dirs = <*_outputs>;

foreach my $in_dir (@dirs) {
	print $in_dir, "\n";
	#system "cd $in_dir";
	my $abs_in_dir = $dir . '/' . $in_dir;
	chdir $abs_in_dir;
	system "cp *_1_GEM_events.bed $abs_resultdir";
}
