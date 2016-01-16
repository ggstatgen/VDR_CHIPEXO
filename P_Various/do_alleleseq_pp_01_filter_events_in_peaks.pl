#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#get only instances in peaks

#inputs
#1 bed with encode blacklist regions (b37?)
#interestinghets.txt file

#file structured as follows
#chrm    snppos          ref     mat_gtyp        pat_gtyp        c_gtyp  phase   mat_all pat_all cA      cC      cG      cT      winning SymCls  SymPval BindingSite     cnv
#1       797440  T       C       T       Y       HETERO  None    None    0       0       0       8       ?       Asym    0.0078125       0       1.0
#1       1314015 C       Y       C       Y       HETERO  None    None    0       0       0       8       ?       Asym    0.0078125       0       1.0
#1       1371459 A       R       G       R       HETERO  None    None    0       0       6       0       ?       Asym    0.03125 0       1.0
#1       1509156 A       A       G       R       HETERO  None    None    0       0       5       0       ?       Sym     0.0625  0       1.0
#1       1706136 T       Y       T       Y       HETERO  None    None    0       0       0       6       ?       Asym    0.03125 0       1.0

my $infile;
GetOptions(
        'i=s'        =>\$infile,
);
if(!$infile){
	print "USAGE: do_alleleseq_pp_01_filter_events_in_peaks.pl -i=<INFILE>\n";
	print "<INFILE> file interestingHets.txt from AlleleSeq\n";
    exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	if($_ =~ /^chrm/){
		print $_, "\n";
		next;
	}
	my $bindingSite = (split /\t/)[16];
	next if($bindingSite ne '1');
	print $_, "\n";
}
close $instream;