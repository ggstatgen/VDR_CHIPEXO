#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#28/4/2014
#Use this after IMPUTE2 and do_impute2_gen_to_vcf_cluster.sh
#You will have had a vcf file with NO sample names and NO chromosome names in part
#eg
##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample_1        sample_2        sample_3        sample_4        sample_5        sample_6        sample_7        sample_8        sample_9        sample_10       sampl
#NA      13302   rs180734498     C       T       .       .       .       GT:GP   ./.:0.605,0.39,0.004    ./.:0.888,0.111,0       #0/0:0.938,0.061,0       ./.:0.835,0.163,0.002   ./.:0.856,0.144,0       0/0:0.911,0.087,0.003   0/0:0.909,0.0

#This file will replace NA with chrX and sample1, etc with a list of samples (taken from file)

my $infile; 
my $samplelist;
GetOptions(
	'i=s'        => \$infile,
	'sample=s'   => \$samplelist
);
if(!$infile){
     print "USAGE: do_impute2_finalize_out_vcf.pl -i=<INFILE>-sample=<SAMPLE_FILE>\n";
     print "<INFILE> .vcf file to process\n";
     print "<SAMPLE_FILE> file with sample names, one per row\n";
     print "ATTENTION: the script assumes the chr name is in <INFILE>\n";
     exit 1;
}
if(!$samplelist){
     print "USAGE: do_impute2_finalize_out_vcf.pl -i=<INFILE> -sample=<SAMPLE_FILE>\n";
     print "<INFILE> .vcf file to process\n";
     print "<SAMPLE_FILE> file with sample names, one per row\n";
     print "ATTENTION: the script assumes the chr name is in <INFILE>\n";
     exit 1;
}
my ($filename, $directory) = fileparse($infile);
my $chr;
if( ($filename =~ /(chr\d+)\.vcf/) || ($filename =~ /(chr[X|Y|M])\.vcf/) ){
	$chr = $1;
}else{
	print "ERROR: unable to identify the chromosome number string from the file name. Aborting.\n";
	exit -1;
}
$filename =~ s/(.*)\..*/$1/;
my $outfile =  $directory . "\/" . $filename . "_FINAL.vcf";

#slurp sample names
my %samplenames;
open (my $instream,  q{<}, $samplelist) or die("Unable to open $samplelist : $!");
while(<$instream>){
	chomp;
	$samplenames{$_} = 1;
}
close $instream;

#build sample string
my @samplestring;
foreach my $sample (sort keys %samplenames){
	push(@samplestring, $sample);
}
my $samplestring = join("\t", @samplestring);

open ($instream,  q{<}, $infile) or die("Unable to open $infile : $!");
open (my $outstream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
while(<$instream>){
	#meta-info / header
	if($_ =~ /^\#\#/){
		print $outstream $_;
		next;
	}
	if($_ =~ /^\#CHROM/){
		my $new_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" . $samplestring;
		print $outstream $new_header, "\n";
		next;
	}
	if($_ =~ /^NA/){
		$_ =~ s/^NA/$chr/;
		print $outstream $_;
		next;
	}else{
		print "This line is not recognised: $_\n. Aborting.";
		exit -1;
	}
	
}
close $instream;
close $outstream;