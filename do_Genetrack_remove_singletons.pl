#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#27/02/2013
#I need to remove all lines of Genetrack output where singletons are found
#a singleton is indicated by genetrack with a std dev = 0.0
#eg 
#chr1    genetrack       .       568373  568383  6.24515535516   +       .       readcount=6;ID=P568378;stddev=2.35702260396;height=6.24515535516
#chr1    genetrack       .       568477  568487  6.0292452985    +       .       readcount=5;ID=P568482;stddev=0.0;height=6.0292452985
#input genetrack gff file
#output genetrack gff file without singletons

my $infile;
GetOptions(
	'input=s'  =>\$infile,
);

if(!$infile){
     print "USAGE: do_Genetrack_remove_singletons.pl -input=<INFILE>\n";
     print "<INFILE> input fastq file\n";
     exit 1;
}
my $outfile = $infile;
$outfile =~ s/(.*)\..*/$1/;
$outfile = $outfile . "_nosgtons.gff";
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
open (my $outstream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
while(<$instream>){
	print $outstream $_ unless($_ =~ /stddev=0.0;/);
}
close $instream;
close $outstream;
