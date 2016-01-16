#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $input;
GetOptions(
        'i=s'		   =>\$input,
);

if(!$input){
     print "USAGE: Brigida.pl -i=<FILE_INPUT>\n";
     print "<FILE INPUT> full path to the input file";
     exit 1;
}

my %data_struct;
open (my $instream,     q{<}, $input) or die("Unable to open $input : $!");
while(<$instream>){
	chomp;
	my @fields = split(" ", $_);
	my @identifiers =  split("_", $fields[0]);

	my $ceppo_ID = shift @identifiers;
	my $sec_ID = shift @identifiers;
	
	#assuming "midasin" or "Dynein_heavy_chain,_N-terminal_region_1" is a signature for the group of ten
	my $signature = join("_", @identifiers);
	$data_struct{$signature}{$ceppo_ID} += 1;
}
close $instream;

foreach my $signature (sort keys %data_struct){
	my %dict_ceppi = (
		"X4000"  => 0, 
		"X4002"  => 0,
		"X4009"  => 0,
		"X4025"  => 0,
		"X4028"  => 0,
		"X4031"  => 0,
		"X4033"  => 0,
		"X4035"  => 0,
		"X4037"  => 0,
		"X4040"  => 0
		);
	foreach my $ceppo (sort keys %{ $data_struct{$signature}} ){
		if(undef $dict_ceppi{$ceppo}){
			print "WARNING: the ceppi dictionary does not contain the ceppo $ceppo.\n";
			print "Are you sure the dictionary contains all the ceppi? ABORTING.\n";
			exit -1;
		}
		$dict_ceppi{$ceppo} = 1;
	}
	foreach my $item (sort keys %dict_ceppi){
		if($dict_ceppi{$item} == 0){
			print "Group: $signature. Ceppo $item is not observed in this group.\n";
		}
	}
}
