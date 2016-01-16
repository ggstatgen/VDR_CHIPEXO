#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#converts a npk file to bed
#if necessary, sorts by pval, cuts the best n, and returns an .npk and a .bed with only those


#10/04
#Testing SPP on chip-exo files. Control is a dummy control with 100 reads selected at random from GM06986
#I want to see if the peaks called by SPP are comparable with others. But IGV won't accept SPP output peaks (narrowpeak)
#I convert them to bed
#23/04 also needed by all the other encode scripts

#from
#chr1    121485104       121485268       .       0       .       301.84984942945 -1      3.11182238802035        82
#to
#chr1    121485104       121485268       x

my $infile;
my $nbest;
GetOptions(
	'i=s'  =>\$infile,
	'b=i'  =>\$nbest
);
if(!$infile){
     print "USAGE: narrowPeak_to_bed.pl -i=<FILE> -b=<N_BEST>\n";
     print "<FILE> input narrowPeak file\n";
     print "<N_BEST> number of best peaks (by p-value) to return (optional, default = all)\n";
     exit 1;
}

if($nbest){
	my $peak_num = `cat $infile | wc -l`;
	if($nbest > $peak_num){
		print STDERR "WARNING:  $infile contains less peaks than you indicated via -b. Ignoring -b=$nbest\n";
	}else{
		#create a new .narrowPeak file
		my($basename, $directory) = fileparse($infile);
		$basename =~ s/(.*)\..*/$1/;
		my $npk_out = $directory . $basename . '_best' . $nbest . '.narrowPeak';
		system "sort -k8nr $infile | head -n $nbest | sort -k1,1V -k2,2g > $npk_out";
		$infile = $npk_out;		
	}

}

open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	my ($chrom, $chromStart, $chromEnd,$name,$score) = (split /\t/)[0, 1, 2,3,4];
	my $newline = $chrom . "\t" . $chromStart . "\t" . $chromEnd . "\t" . $name . "\t" . $score; 
	print $newline, "\n";
}
close $instream;
