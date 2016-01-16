#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#produce a bed file of the funseq output. You can produce different ones based on
#all variants
#recur
#vdr peak
#recur + vdr peak

my $infile;
my $FILTER_VDR_ONLY;
my $FILTER_RECUR;
GetOptions(
        'i=s'        =>\$infile,
        'vdr'   	 =>\$FILTER_VDR_ONLY,
        'recur'   	 =>\$FILTER_RECUR,
);

if(!$infile){
	print "USAGE: do_funseq_output2bed.GAT.pl -i=<INFILE> -vdr -recur\n";
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
	next if($_ =~ /^\#/);
	next if($_ eq '');

	#get info field
	my ($chr, $pos, $rsID, $info) = (split /\t/)[0,1,2,7];
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
	my $id = $chr . '-' . $pos . '-' . $rsID; #recurrent is  printed once per sample in output.vcf
	$id_to_info{$id} = $info;
}

#print bed
foreach my $item (sort keys %id_to_info){
	my @item = split('-', $item);
	print $item[0] . "\t" . ( $item[1] - 1) . "\t" . $item[1] . "\n";
}