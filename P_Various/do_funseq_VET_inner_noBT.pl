#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;
use List::Util qw(sum);

#version of do_funseq_ET_02_freqmatch_bootstrap_multithreading_inner.pl which does NOT use bedtools for the fg/bg intersection. Faster?

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
#$input_bg = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/ET-VDRBV-REP/LDCLUMP_REP_FG_VDRBVs_Output_noDBRECUR_TEMP_bgmatrix.txt";
#$job_counter = 3;
#$ph_map = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/ET-VDRBV-REP/LDCLUMP_REP_FG_VDRBVs_Output_noDBRECUR_TEMP_diseaseblocks_phmap.txt";
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

my %test_vars;
my @METRIC_BG;
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

###################
#2 get phenotype interval data structure
###################
my %ph_intervals; #hash of arrays
open ($temp_stream,  q{<}, $out_temp_buffer_ph_bed_processed) or die("Unable to open $out_temp_buffer_ph_bed_processed : $!");
while(<$temp_stream>){
	chomp;
	my ($chr, $start, $stop) = split("\t",$_);
	my $interval = $start . '-' . $stop;
	push(@{$ph_intervals{$chr}}, $interval);
}
close $temp_stream;

####################
#3 get foreground data & metrics
#####################
tie *FILE,   'IO::Zlib', $input_fg, "rb";
while(<FILE>){
	chomp;
	my ($chr,$pos) = (split /\t/)[0,2];
	push(@{$test_vars{$chr}}, $pos);
}
close FILE;

#see if the chr of the fg variant is in the gwas set.
#if no move on to next chr;
#if yes, get position of the fg variant and check against all gwas intervals to see if it's in any
my $METRIC_FG = 0;
foreach my $chr (sort keys %test_vars){
		if($ph_intervals{$chr}){
			foreach my $pos (@{$test_vars{$chr}}) {
	    		foreach my $interval (@{$ph_intervals{$chr}}) {
	        		my ($start, $stop) = split('-',$interval);
	        		if($pos >= $start && $pos <= $stop){
	        			$METRIC_FG++;
	        			last;# each variant can be in only one interval (intervals are unique)
	        		}
	    		}
	    	}
		}		
}
if($METRIC_FG eq '0'){
	print STDERR "ERROR: for phenotype: $ph_name - by definition this number should not be zero: $METRIC_FG. Aborting..\n";
	exit -1;
}

###################
#4 slurp background matrix into array of arrays
###################

open ($temp_stream,  q{<}, $input_bg) or die("Unable to open $input_bg : $!");
while(<$temp_stream>){
	chomp;
	my @row = split("\t",$_);
	
	my $METRIC_BG = 0;
	foreach my $variant (@row){
		my ($chr, $pos) = split('-', $variant);	
		if($ph_intervals{$chr}){
    		foreach my $interval (@{$ph_intervals{$chr}}) {
        		my ($start, $stop) = split('-',$interval);
        		if($pos >= $start && $pos <= $stop){
        			$METRIC_BG++;
        			last;# each variant can be in only one interval (intervals are unique)
        		}
    		}
		}					
	}
	push(@METRIC_BG, $METRIC_BG);
}
close $temp_stream;

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