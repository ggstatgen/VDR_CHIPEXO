#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#just like do_funseq_GEUVADIS_to_bed.pl, but name will only be [gene|exon|trratio]
#if you want to do a GAT separately by type (exon, gene, ttratio)
#results are in
#/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/d_GAT/GAT_eQTL_geuvadis_YRI_gene_exon_ttratio_all_bytype_r10000
#significance

my $infile;
my $pop; #ceu/yri
my $type; #geneQTL/exonQTL/trratio
GetOptions(
        'i=s'      =>\$infile,
        'pop=s'    =>\$pop,
        'type=s'   =>\$type
);
if(!$infile){
	print "USAGE: do_funseq_GEUVADIS_to_bed.pl -i=<INFILE> \n";
    print "<INFILE> GEUVADIS GZIPPED tsv file\n";
    
    exit 1;
}
my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
my $output  = $directory . $basename . '.bed';

if($basename =~ /gene/){
	$type = 'geneQTL';
}elsif($basename =~ /exon/){
	$type = 'exonQTL';
}elsif($basename =~ /trratio/){
	$type = 'trratioQTL';
}else{
	print STDERR "Error: impossible to determine QTL type from filename: $basename. Aborting.\n";
	exit -1;
}

open (my $outstream,  q{>}, $output) or die("Unable to open $output : $!");
tie *FILE,   'IO::Zlib', $infile, "rb";
while (<FILE>)	{ 
	chomp;
	next if($_ =~ /^SNP_ID/);
	my @fields = split("\t", $_);
	
	#skip indels
	if($fields[6] =~ /\./){
		print STDERR "Warning: this appears to be an indel: $fields[6]. Skipping..\n";
		next;
	}
	print $outstream $fields[4] . "\t" . ($fields[6]-1) . "\t" . $fields[6] . "\t" . $type . "\n"; 
}
close FILE;
close $outstream;

