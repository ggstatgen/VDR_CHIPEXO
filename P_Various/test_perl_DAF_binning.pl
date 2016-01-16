#!/usr/bin/perl
use strict;
use warnings;
use Statistics::R;

my $BIN_NUMBER = 20;

my @DAF_vals = (0.266845129, 0.353675275,0.030585143,0.791285527,0.60805832,0.433466716,0.039369439,0.48918122,0.740626586,0.2,0.21,0.22);

my %DAF_binning;
foreach my $item (@DAF_vals){
	my $bin_id = int($item * $BIN_NUMBER ) / $BIN_NUMBER;
	#$DAF_binning{$bin_id}{'test'}++;
	$DAF_binning{$bin_id}++;
}

for (0..20){
    $_ = $_/20;
    
    #if(ref $DAF_binning{$_}){
     if($DAF_binning{$_}){  	
    	#get number of elements in this bin
    	#my $number_in_this_bin = scalar(keys %{$DAF_binning{$_}});
    	my $number_in_this_bin = $DAF_binning{$_};
    	print "$_\t" .$number_in_this_bin . "\n";
    }else{
    	print "$_\t0\n";
    }	
}

