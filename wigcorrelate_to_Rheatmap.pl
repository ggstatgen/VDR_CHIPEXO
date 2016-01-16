#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#I want to create an input suitable for R (heatmap)
#This program accepts one text file in the format
#<PATH>/GM06986_logLR.bw  <PATH>/GM06989_logLR.bw  0.603199
#<PATH>/GM06986_logLR.bw  <PATH>/GM06997_logLR.bw  0.372119

#data in the input file is as follows:
#given n datapoints A B C D, you have:

#A B x_1
#A C x_2
#A D x_3
#B C x_4
#B D x_5
#C D x_6
#for each line split 1st, 2nd and rho
#create hash for all unique datasets. Key = dataset. Value = 1. %data = {A->1, B->1, ..., Z->1}
#create ordered hash for data. Key = sort(pair_1, pair_2), Value = rho
#include all self data with value: Key = (pair_i,pair_i), value = 1

my $infile; 
GetOptions(
	'i=s'  => \$infile,
);
if(!$infile){
     print "USAGE: wigcorrelate_to_R_heatmap.pl -i=<INFILE>\n";
     print "<INFILE> text file obtained with wigcorrelate\n";
     exit 1;
}
my ($filename, $directory) = fileparse($infile);
$filename =~ s/(.*)\..*/$1/;
my $outfile =  $directory . "\/" . $filename . "_R.dat";

open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
my %data_label_hash;
my %data_points_hash;
while(<$instream>){
	my $filename_i; my $filename_j; my $directory;
	my ($data_i, $data_j, $rho) = split(" ", $_);
	
	#processing on the strings to remove path, etc
	($filename_i, $directory) = fileparse($data_i);
	($filename_j, $directory) = fileparse($data_j);
	$filename_i =~ s/(.*)\_.*/$1/;
	$filename_j =~ s/(.*)\_.*/$1/;
	
	#collect sample names
	$data_label_hash{$filename_i} = 1;
	$data_label_hash{$filename_j} = 1;
	#collect data
	$data_points_hash{join("_", sort($filename_i, $filename_j))} = $rho;
	$data_points_hash{join("_", $filename_i,$filename_i)} = 1;
	$data_points_hash{join("_", $filename_j,$filename_j)} = 1;
}
close $instream;

open (my $outstream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
#print header
print $outstream " ";
foreach my $item_i (sort keys %data_label_hash){
	print $outstream $item_i, " ";
}
print $outstream "\n";

my $dataline;
my @dataline;
foreach my $item_i (sort keys %data_label_hash){
	push(@dataline, $item_i);
	foreach my $item_j (sort keys %data_label_hash){
		if(defined $data_points_hash{join("_", sort($item_i,$item_j))}){
			push(@dataline, $data_points_hash{join("_", sort($item_i,$item_j))});
		}else{
			push(@dataline, 'NA');
		}		
	}
	print $outstream join(" ", @dataline), "\n";
	@dataline = ();
}
close $outstream;