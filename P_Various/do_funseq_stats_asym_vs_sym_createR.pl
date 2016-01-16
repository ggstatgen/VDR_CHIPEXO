#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#statistical test between VDR-BV and SYM variants.
#you NEED to know whether SYM are LOB or GOB
#in the figure, you assigned LOB or GOB to SYM in R, after computer alt+1/ref+1

#ultimately you want to do SYM vs ASYM based on the magnitude of the effect

#input format
#type    position        ref     alt
#MOTIFBR pos_1   1       11
#LOB
#GOB
#MOTIFBR_SYM
#LOB_SYM
#GOB_SYM

my $infile;
GetOptions(
        'i=s'	=>\$infile
);
#$infile = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/vdrrxr_phenotypes/d_FIGURE3_STATS/asym_vs_sym.Rdata";

if(!$infile){
	print "USAGE: do_funseq_stats_asym_vs_sym_createR.pl -i=<INFILE>\n";
	print "<INFILE> .Rdata file with phenotypes and positions for sym and asym";
    exit 1;
}
print STDOUT "type\tposition\tref\talt\tmetric\n";
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	next if($_ =~ /^type/);
	
	chomp;
	my $metric;
	my ($type, $pos, $ref, $alt) = split("\t", $_);
	#if($type eq 'MOTIFBR'){ $type = 'LOB'; }
	
	if($type =~ 'SYM'){
		if($type =~ 'MOTIFBR'){ $type = 'LOB'; }
		my $label = 'SYM';
		
		if( ($ref eq 'NA')  && ($alt eq 'NA') ){
			print STDOUT $label . "\t" . $pos . "\t" . $ref . "\t" . $alt . "\t" .  'NA' . "\n";
			next;
		}	
		
		if($type =~ 'LOB'){
			$metric = ( $ref + 1 ) / ( $alt  + 1 );
			print STDOUT $label . "\t" . $pos . "\t" . $ref . "\t" . $alt . "\t" .  $metric . "\n";
		}elsif($type =~ 'GOB'){
			$metric = ( $alt + 1 ) / ( $ref  + 1 );
			print STDOUT $label . "\t" . $pos . "\t" . $ref . "\t" . $alt . "\t" .  $metric . "\n";		
		}else{
			print STDERR "Error: type: $type not recognised. Aborting..\n";
			exit -1;
		}			
	}else{
		if($type eq 'MOTIFBR'){ $type = 'LOB'; }
		my $label = 'ASYM';
		
		if( ($ref eq 'NA')  && ($alt eq 'NA') ){
			print STDOUT $label . "\t" . $pos . "\t" . $ref . "\t" . $alt . "\t" .  'NA' . "\n";
			next;
		}	
		
		if($type eq 'LOB'){
			$metric = ( $ref + 1 ) / ( $alt  + 1 );
			print STDOUT $label . "\t" . $pos . "\t" . $ref . "\t" . $alt . "\t" .  $metric . "\n";
		}elsif($type eq 'GOB'){
			$metric = ( $alt + 1 ) / ( $ref  + 1 );
			print STDOUT $label . "\t" . $pos . "\t" . $ref . "\t" . $alt . "\t" .  $metric . "\n";		
		}else{
			print STDERR "Error: type: $type not recognised. Aborting..\n";
			exit -1;
		}		
	}
}
close $instream;
