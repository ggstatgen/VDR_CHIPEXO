#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#subset the output vcf coming out of funseq in subsets
#-vdr peaks
#-recurs
#-vdr peaks and recurs

my $infile;
my $FILTER_VDR_ONLY;
my $FILTER_RECUR;
GetOptions(
        'i=s'        =>\$infile,
        'vdr'   	 =>\$FILTER_VDR_ONLY,
        'recur'   	 =>\$FILTER_RECUR,
);
if(!$infile){
	print "USAGE: do_funseq_output2vcfsubset.pl -i=<INFILE> -vdr -recur\n";
    print "<INFILE> vcf output of funseq2\n";
    print "(optional)<vdr> whether to only use ASB variants in a VDR peak (default=no)\n";
    print "(optional)<recur> whether to only use ASB variants which recur (default=no)\n";
    exit 1;
}

#get all snps or the subset with VDR PEAK, RECUR or BOTH annotation
my %id_to_info;
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	if($_ =~ /^\#/){ print $_, "\n"; next; }
	next if($_ eq '');

	#get info field
	my @all = split("\t", $_);
	my $info = $all[7];
	#filters-------------------------------------------------
	if($FILTER_RECUR){
		next unless($info =~ /RECUR/);
	}
	if($FILTER_VDR_ONLY){
		my $info_ncenc;
		my @info = split(";", $info);
		foreach my $item (@info){
			$info_ncenc = $item if($item =~ /^NCENC/);
		}
		next if(!$info_ncenc);
		next unless($info_ncenc =~ /TFP(.*)VDR/);
	}
	#-------------------------------------------------------
	$all[0] =~ s/chr(.*)/$1/; 
	print $all[0] . "\t" . $all[1] . "\t" . $all[2] . "\t" . $all[3] . "\t" . $all[4] . "\t" . $all[5] . "\t" . $all[6] . "\tSEE_FUNSEQ_MASTER_FILE" . "\n";
	#print $_, "\n";
}

