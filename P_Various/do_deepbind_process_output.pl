#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#process the output of deepbind so that the header reports the TF name, not the code.
#get the code from the table produced by deepbind with the command
#deepbind --dump-info example.ids

#sample out



my $in_deepbind;
my $in_info;

GetOptions(
        'out=s'       =>\$in_deepbind,
        'info=s'      =>\$in_info
);

#$in_deepbind = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/DEEPBIND/results_cI.txt";
#$in_info = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/DEEPBIND/info.tsv";

if(!$in_deepbind){
	print "USAGE: do_deepbind_process_output.pl -out=<DEEPBIND_OUT> -info=<DEEPBIND_INFO>\n";
    print "<DEEPBIND_OUT> tsv file of scores produced by deepbind\n";
    print "<DEEPBIND_INFO> info map of ids to tfs produced by deepbind\n";
    exit 1;
}
if(!$in_info){
	print "USAGE: do_deepbind_process_output.pl -out=<DEEPBIND_OUT> -info=<DEEPBIND_INFO>\n";
    print "<DEEPBIND_OUT> tsv file of scores produced by deepbind\n";
    print "<DEEPBIND_INFO> info map of ids to tfs produced by deepbind\n";
    exit 1;
}

my($basename, $directory) = fileparse($in_deepbind);
$basename =~ s/(.*)\..*/$1/;
my $output = $directory . $basename . '.Rdata';

#get id>tf map
#ID	Protein	Type	Species	Family	Class	Experiment	Experiment Details	Model	Cite	Labels	Path	Comment
#D00210.001	RBFOX1	RBP	Homo sapiens	RRM		RNAcompete	['RNAcompeteID=RNCMPT00168']	deepbind 0.1	PMID 23846655		E:/Results/deepbind/rnac/AB/final/RNCMPT00168	

my %id2tf;
open (my $instream,   q{<}, $in_info) or die("Unable to open $in_info : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^ID/);
	my ($id, $tf) = (split /\t/)[0,1];
	$id2tf{$id} = $tf;
}
close $instream;

#create new map 

open ($instream,   q{<}, $in_deepbind) or die("Unable to open $in_deepbind : $!");
open (my $outstream,   q{>}, $output) or die("Unable to open $output : $!");
while(<$instream>){
	chomp;
	
	#header
	if($_ =~ /^D/){ #the IDs start with D it seems
		my @header = split("\t", $_);
		my @processed_header = map { $id2tf{$_} } @header;
		print $outstream join("\t",@processed_header) . "\n";
		next;
	}
	print $outstream $_, "\n";
}
close $instream;
close $outstream;
