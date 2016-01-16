#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#downloaded from here and adapted
#http://info.gersteinlab.org/ACT_Tool

my $infile;
GetOptions(
        'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_sgl_to_wig.pl -i=<FILE>\n";
     print "<FILE> input sgl file\n";
     exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");

my $chr="chr";
my $start=-1;
my $sig=0;
while(<$instream>)
{
	my @line = split/\s+/, $_;
	my $temp_chr = $line[0];
	my $temp_start = $line[1]-1;
	my $temp_sig = $line[2];	
	if ($chr ne $temp_chr){
		$chr = $temp_chr;
		$start = $temp_start + 1;
		$sig = $temp_sig;
	}
	else{
		print "$chr\t$start\t$temp_start\t$sig\n";
		$start = $temp_start + 1;
		$sig = $temp_sig;
	}
	
}
close $instream;
