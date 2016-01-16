#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;
use List::Util qw(sum);

#obsolete, slower than the other one where I do not use bedtools

#trying to parallelize the most time consuming part of do_funseq_ET_02_freqmatch_bootstrap.pl: the enrichment calculations
#I run this perl script, 1 time per phenotype, from within do_funseq_ET_02_freqmatch_bootstrap.pl
#the outcome is also collected by do_funseq_ET_02_freqmatch_bootstrap.pl, one line per instance, for the calculation of the corrected pvalue.

#inputs
#phenotype name string
#file with foreground vars
#file with matrix of 1000 x sets of background vars

my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";

my $input_fg;
my $input_bg;
my $ph_map;
my $blocks_gwas;
my $job_counter;
my $ph_name;
GetOptions(
        'fg=s'      =>\$input_fg,
        'bg=s'      =>\$input_bg,
        'ph=s'      =>\$ph_map,
        'n=i'       =>\$job_counter,
        'gwas=s'    =>\$blocks_gwas
);

#$input_fg = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/ET-VDRBV-REP/LDCLUMP_REP_FG_VDRBVs_Output_noDBRECUR.bed.gz";
#$input_bg = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/ET-VDRBV-REP/#LDCLUMP_REP_FG_VDRBVs_Output_noDBRECUR_TEMP_bgmatrix.txt";
#$ph_name = '2 hour glucose';
#$blocks_gwas = "/net/isi-scratch/giuseppe/indexes/GWAS_GRASP/LD_PLINK_CEU_GRASP2_plus_Beecham2013.bed.gz";

if(!$input_fg){
	print STDERR "USAGE: [script.pl] -fg=<FGfile> -bg=<BGmatrix> -ph=<phmap> -n=<JOB_ID>  -ld=<GWAS_LD>\n";
    exit 1;
}
if(!$input_bg){
	print STDERR "USAGE: [script.pl] -fg=<FGfile> -bg=<BGmatrix> -ph=<phmap> -n=<JOB_ID>  -ld=<GWAS_LD>\n";
    exit 1;
}
if(!$ph_map){
	print STDERR "USAGE: [script.pl] -fg=<FGfile> -bg=<BGmatrix> -ph=<phmap> -n=<JOB_ID>  -ld=<GWAS_LD>\n";
    exit 1;
}
if(!$job_counter){
	print STDERR "USAGE: [script.pl] -fg=<FGfile> -bg=<BGmatrix> -ph=<phmap> -n=<JOB_ID>  -ld=<GWAS_LD>\n";
    exit 1;	
}
if(!$blocks_gwas){
	print "USAGE: [script.pl] -fg=<FGfile> -bg=<BGmatrix> -ph=<phmap> -n=<JOB_ID>  -ld=<GWAS_LD>\n";
    exit 1;
}
my($basename, $directory) = fileparse($input_fg);
$basename =~ s/(.*)\..*/$1/;
$basename =~ s/(.*)\..*/$1/;

my $out_temp_buffer_ph_bed           = $directory . 'TEMP_phenotype_'         .  $job_counter . '.bed';
my $out_temp_buffer_ph_bed_processed = $directory . 'TEMP_phenotype_'         .  $job_counter . 'processed.bed';
my $out_temp_input_bg                = $directory . 'TEMP_phenotype_'         .  $job_counter . '_bgset.bed';

my @METRIC_BG = ();

#############
#0 get phenotype name from map
#############
open (my $temp_stream,  q{<}, $ph_map) or die("Unable to open $ph_map : $!");
while(<$temp_stream>){
	chomp;
	my ($id, $this_ph_name) = split("\t", $_);
	if($id eq $job_counter){
		$ph_name = $this_ph_name;
	}
}
close $temp_stream;

unless($ph_name){
	print STDERR "phenotype ID: $job_counter - could not identify a phenotype name from the map.Aborting\n";
	exit -1;
}

#############
#1 get bed of GWAS LD blocks for the phenotype in $ph_name only
#############
open ($temp_stream,  q{>}, $out_temp_buffer_ph_bed) or die("Unable to open $out_temp_buffer_ph_bed : $!");
tie *FILE,   'IO::Zlib', $blocks_gwas, "rb";
while(<FILE>){
	chomp;
	next if($_ eq '');
	my ($ld_block_chr, $ld_block_start, $ld_block_stop, $phenotype) = (split /\t/)[0,1,2,11];

	if( lc($phenotype) eq lc($ph_name) ){
		my $line = $ld_block_chr . "\t" . $ld_block_start . "\t" . $ld_block_stop  . "\t" . $ph_name;
		print $temp_stream $line, "\n";
	}
}
close FILE;
close $temp_stream;

system "sort -k1,1V -k2,2g $out_temp_buffer_ph_bed | $BEDTOOLS merge -i stdin > $out_temp_buffer_ph_bed_processed";

####################
#2 get foreground metrics
#####################
my $METRIC_FG = `$BEDTOOLS intersect -a $input_fg -b $out_temp_buffer_ph_bed_processed | wc -l`;
chomp $METRIC_FG;
if($METRIC_FG eq '0'){
	print STDERR "ERROR: for phenotype: $ph_name - by definition this number should not be zero: $METRIC_FG. Aborting..\n";
	exit -1;
}
###################
#3 slurp background matrix into array of arrays
###################
my @random_sample_matrix;
open ($temp_stream,  q{<}, $input_bg) or die("Unable to open $input_bg : $!");
while(<$temp_stream>){
	chomp;
	 push @random_sample_matrix, [ split /\t/ ];
}
close $temp_stream;

#############
#4 generate background beds, intersect each with this phenotype's bed and save all intersections
#############
foreach my $row (@random_sample_matrix) { #each row is a randomized sample of #foreground SNPs
	open (my $tempfile,  q{>}, $out_temp_input_bg) or die("Unable to open $out_temp_input_bg : $!");
	foreach my $col (@{$row}) {
		my ($chr, $pos) = split('-', $col);
  		print $tempfile $chr . "\t" . ($pos - 1) . "\t" . $pos . "\n";
	}
	close $tempfile;
	#intersect
	my $METRIC_BG = `sort -k1,1V -k2,2g $out_temp_input_bg | $BEDTOOLS intersect -a stdin -b $out_temp_buffer_ph_bed_processed | wc -l`;
	unlink $out_temp_input_bg;
	chomp $METRIC_BG;
	push(@METRIC_BG, $METRIC_BG);
}
##############
#5 stats
##############
#now you have, for this phenotype, an intersection for the foreground and 1000 intersections for the background
my ($obs, $exp, $fc, $l2fc, $pval) = get_stats($METRIC_FG,\@METRIC_BG);
my $stats = $obs . "\t" . $exp . "\t" . $fc . "\t" . $l2fc . "\t" . $pval;

print STDOUT  $ph_name . "\t" . $stats; #this needs to be captured

unlink $out_temp_buffer_ph_bed;
unlink $out_temp_buffer_ph_bed_processed;


















################
#subs
################
sub get_stats{
	my ($value_fg, $values_bg) = @_;
	my $counter = 0;
	
	foreach my $item (@$values_bg){
		if($item >= $value_fg){
			$counter++;
		}
	}
	my $p =  ( ( 1 + $counter) / (@$values_bg + 1)    );
	#what is the expected value? the mean of the values in the randomisation?
	my $obs = $value_fg;
	my $exp	= mean(@$values_bg);
	my $fold_change;
	#if($exp > 0){
		$fold_change = ( 1.0 + $obs) / ( 1.0 + $exp);
	#}else{
	#	$fold_change = 'inf';
	#}
	return ($obs, $exp, $fold_change, log2($fold_change), $p);
}

sub log2 {
	my $n = shift;
	return log($n)/log(2);
}

sub mean {
    return sum(@_)/@_;
}