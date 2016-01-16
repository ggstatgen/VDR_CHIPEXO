#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#use this to remove all peaks in blacklisted regions from a hg19 bed file
#https://sites.google.com/site/anshulkundaje/projects/blacklists
#you can choose between hg19 and mm9 (and at a push mm10, which I obtained from mm9 via liftover - however this had non-mapping locations)

#this uses bedtools subtract
#Usage:   bedtools subtract [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>
#
#Options: 
#	-f	Minimum overlap required as a fraction of A.
#		- Default is 1E-9 (i.e., 1bp).
#		- (FLOAT) (e.g. 0.50)
#
#	-s	Require same strandedness.  That is, only subtract hits in B
#		that overlap A on the _same_ strand.
#		- By default, overlaps are subtracted without respect to strand.
#
#	-S	Force strandedness.  That is, only subtract hits in B that
#		overlap A on the _opposite_ strand.
#		- By default, overlaps are subtracted without respect to strand.
#
#	-A	Remove entire feature if any overlap.  That is, by default,
#		only subtract the portion of A that overlaps B. Here, if
#		any overlap is found (or -f amount), the entire feature is removed.

my $ENCODE_DIR = '/net/isi-scratch/giuseppe/indexes/ENCODE_BLACKLIST';
my $BEDTOOLS = '/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin';
my $ENCODE_FILE;

my $dir;my $ref; my $ext; my $outex;
my $resultdir = 'd_bed_RBL';
my @files;
GetOptions(
	'd=s'    =>\$dir,
	'r=s'    =>\$ref,
	'ext=s'  =>\$ext,
);
if(!$dir){
     print "USAGE:  do_bed_remove_blacklisted_peaks.pl -d=<DIRECTORY> -r=<REFERENCE> -ext=<EXTENSION>\n";
     print "<DIRECTORY> where the bed files are\n";
     print "<REFERENCE> genome to use [hg19|b37|mm9|mm10]\n";
     print "<EXTENSION> [b|n|g] (.bed,. narrowPeak, .gappedPeak - default: bed)\n";
     exit 1;
}
if(!$ref){
     print "USAGE:  do_bed_remove_blacklisted_peaks.pl -d=<DIRECTORY> -r=<REFERENCE> -ext=<E>\n";
     print "<DIRECTORY> where the bed files are\n";
     print "<REFERENCE> genome to use [hg19|b37|mm9|mm10]\n";
     print "<EXTENSION> [b|n] (for bed or narrowPeak - default: bed)\n";
     exit 1;
}

#check $ref
if($ref eq 'hg19'){
	$ENCODE_FILE = $ENCODE_DIR . '/hg19-wgEncodeDacMapabilityConsensusExcludable.bed';
}elsif($ref eq 'mm9'){
	$ENCODE_FILE = $ENCODE_DIR . '/mm9-blacklist.bed';
}elsif($ref eq 'mm10'){
	$ENCODE_FILE = $ENCODE_DIR . '/mm10-blacklist_liftover.bed';
}elsif($ref eq 'b37'){
	$ENCODE_FILE = $ENCODE_DIR . '/b37-wgEncodeDacMapabilityConsensusExcludable.bed';
}else{
	print "-r argument not recognised.\n";
	exit -1;
}
if(!$ext){
        @files = <*.bed>;
        $resultdir = 'd_bed_RBL';
	$outex = ".bed";
}elsif($ext eq 'b'){
        @files = <*.bed>;
        $resultdir = 'd_bed_RBL';
	$outex = ".bed";
}elsif($ext eq 'n'){
        @files = <*.narrowPeak>;
        $resultdir = 'd_npk_RBL';
	$outex = ".narrowPeak";
}elsif($ext eq 'g'){
	@files = <*.gappedPeak>;
	$resultdir = 'd_gpd_RBL';
	$outex = ".gappedPeak";
}else{
        print "$ext not recognised. Exiting...\n";
        exit -1;
}

chdir $dir;
system "mkdir $resultdir";
my $abs_resultdir = $dir . '/' . $resultdir;

foreach my $file (@files) {
	my $outfile = $file;
	$outfile =~ s/(.*)\..*/$1/;

	my $infile = $dir . "/" . $file;
	$outfile = $abs_resultdir . "/" .  $outfile .   "_RBL" . $outex;
	#print $infile . "\t" . $outfile . "\n";
	system "$BEDTOOLS/bedtools subtract -A -a $infile -b $ENCODE_FILE > $outfile";
}


#my @dirs = <*_outputs>;
#
#foreach my $in_dir (@dirs) {
#	print $in_dir, "\n";
#	#system "cd $in_dir";
#	my $abs_in_dir = $dir . '/' . $in_dir;
#	chdir $abs_in_dir;
#	system "cp *_2_GEM_events.bed $abs_resultdir";
#}
