#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#statistical test between LOB and GOB
#you need to create a third column with ref+1/alt+1 for lob and alt+1/ref+1 for gob

#input format
#type    position        ref     alt
#MOTIFBR pos_1   1       11
#LOB
#GOB
#SYM

my $infile;
GetOptions(
        'i=s'	=>\$infile
);
if(!$infile){
	print "USAGE: do_funseq_lob_vs_gob_createR.pl -i=<INFILE>\n";
	print "<INFILE> .Rdata file with phenotypes and positions\n";
    exit 1;
}

#my $infile = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/vdrrxr_phenotypes/all.Rdata";

print STDOUT "type\tposition\tref\talt\tmetric\n";
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	next if($_ =~ /^type/);
	
	chomp;
	my $metric;
	my ($type, $pos, $ref, $alt) = split("\t", $_);
	if($type eq 'MOTIFBR'){ $type = 'LOB'; }
	if($type eq 'SYM'){ next; }
	
	if( ($ref eq 'NA')  && ($alt eq 'NA') ){
		print STDOUT $type . "\t" . $pos . "\t" . $ref . "\t" . $alt . "\t" .  'NA' . "\n";
		next;
	}
	
	if($type eq 'LOB'){
		$metric = ( $ref + 1 ) / ( $alt  + 1 );
		print STDOUT $type . "\t" . $pos . "\t" . $ref . "\t" . $alt . "\t" .  $metric . "\n";
	}elsif($type eq 'GOB'){
		$metric = ( $alt + 1 ) / ( $ref  + 1 );
		print STDOUT $type . "\t" . $pos . "\t" . $ref . "\t" . $alt . "\t" .  $metric . "\n";		
	}else{
		print STDERR "Error: type: $type not recognised. Aborting..\n";
		exit -1;
	}
}
close $instream;
