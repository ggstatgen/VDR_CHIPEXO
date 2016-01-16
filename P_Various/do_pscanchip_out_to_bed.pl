#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#pscan chip is a motif finder based on weeder
#http://159.149.160.51/pscan_chip_dev/

#17/07/2014
#add track name on top
#add score in peak name so it gets picked up by UCSC

#9/10/2014
#choose between outputting bed of REGIONS or BINDINGSITES

#I'm using it to find a consensus motif around the het snps I have found to have asymmetric binding with the AlleleSeq pipeline
#It conveniently returns a .ris file containing the intervals with the enriched motifs
#This file is in the following format:
#CHR     REG_START       REG_END REG_STRAND      ABS_SITE_START  ABS_SITE_END    REL_SITE_START  REL_SITE_END    SITE_STRAND     SCORE   SITE
#chr5    131986387       131986536       0       131986469       131986483       7       20      -0      0.965192        TGAACCCTGTGACCT
#chr5    1297635 1297784 0       1297697 1297711 -13     0       -0      0.941996        TGAACTCCATGAACT
#chr3    189207271       189207420       0       189207374       189207388       28      41      0       0.9018  AGGCCATTGAGTTCA
#chr6    158182596       158182745       0       158182674       158182688       3       16      -0      0.894597        TGAACCGTGTGACCT
#chr7    139586182       139586331       0       139586226       139586240       -31     -18     -0      0.890353        GGAACCTATTAACCT

#This scripts extracts a bed file from the ris, for visualisation in igv.
#The format will be
#CHR	 ABS_SITE_START	ABS_SITE_END	SITE(SITE_STRAND)(SCORE)	SCORE
#0,4,5,8,9,10

#you may threshold on the min score wanted

my $infile;
my $MIN_SCORE;
my $interval_type;
GetOptions(
        'i=s'        =>\$infile,
        't=s'        =>\$interval_type,
        'm=f'        =>\$MIN_SCORE
);
if(!$infile){
	print "USAGE: do_pscanchip_out_to_bed.pl -i=<INFILE> -t=<INTERVAL_TYPE> -m=<MIN_SCORE>\n";
    print "<INFILE> file .ris from Pscanchip\n";
    print "<INTERVAL_TYPE> one of [reg|bds] - full peak region or only exact binding site in the bed?\n";
    print "<MIN_SCORE> lower threshold on score (eg 0.8) (default: none)\n";
    print "NOTE: UCSC track name will be derived from file name, so make it INFORMATIVE!\n";
    exit 1;
}
if(!$interval_type){
	print "USAGE: do_pscanchip_out_to_bed.pl -i=<INFILE> -t=<INTERVAL_TYPE> -m=<MIN_SCORE>\n";
    print "<INFILE> file .ris from Pscanchip\n";
    print "<INTERVAL_TYPE> one of [reg|bds] - full peak region or only exact binding site in the bed?\n";
    print "<MIN_SCORE> lower threshold on score (eg 0.8) (default: none)\n";
    print "NOTE: UCSC track name will be derived from file name, so make it INFORMATIVE!\n";
    exit 1;
}
unless($interval_type eq 'bds'  || $interval_type eq 'reg'){
	print "USAGE: do_pscanchip_out_to_bed.pl -i=<INFILE> -t=<INTERVAL_TYPE> -m=<MIN_SCORE>\n";
    print "<INFILE> file .ris from Pscanchip\n";
    print "<INTERVAL_TYPE> one of [reg|bds] - full peak region or only exact binding site in the bed?\n";
    print "<MIN_SCORE> lower threshold on score (eg 0.8) (default: none)\n";
    print "NOTE: UCSC track name will be derived from file name, so make it INFORMATIVE!\n";
    exit 1;
}


#get basename for track name
#eg
#track type=bed name="CREB1" description="pscanchip motifs for CREB1, score > 0.6" db=hg19

my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
my $track_info;
if($MIN_SCORE){
	$track_info = 'track type=bed name="' . $basename . '" description="pscanchip motifs for ' . $basename . ', score > ' . $MIN_SCORE . '"';
}else{
	$track_info = 'track type=bed name="' . $basename . '" description="pscanchip motifs for ' . $basename . '"';
}
print $track_info, "\n";
print STDERR "THRESHOLDING ON SCORE: $MIN_SCORE\n" if($MIN_SCORE);
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
#print track
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^CHR/);
	
	my $chr; my $site_start; my $site_end; my $site_strand; my $score; my $site;
	
	($chr,$site_start,$site_end,$site_strand,$score,$site) = (split /\t/)[0,1,2,8,9,10] if($interval_type eq 'reg');
	($chr,$site_start,$site_end,$site_strand,$score,$site) = (split /\t/)[0,4,5,8,9,10] if($interval_type eq 'bds');
	
	next if (!$chr);
	next if(  $MIN_SCORE && ($score < $MIN_SCORE) );

	my $info = $site . '(' . $site_strand . ')' . '(' . $score . ')';
	print $chr . "\t" . $site_start . "\t" . $site_end . "\t" . $info . "\t" . $score . "\n";	
}
close $instream;
