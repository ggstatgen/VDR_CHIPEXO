#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#I have found interesting hits when using ALL the mapped data, not only the data mapped to peaks
#However I want to remove all interestinghits  outside of peaks which fall in the Encode Blacklist regions

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


#make temp interesting hets file with another column, to turn into a bed. Then, use bedtools. Then remove the column from the output


my $TOOL_PATH       = '/net/isi-scratch/giuseppe/tools';
my $BEDTOOLS        = $TOOL_PATH . '/bedtools-2.22.1/bin';
my $BLACKLIST_REGIONS = "/net/isi-scratch/giuseppe/indexes/BLACKLIST_ALLTRACKS/b37_anshuldac_duke.bed";
my $infile;

GetOptions(
        'i=s'        =>\$infile,
);

#$infile = '/net/isi-scratch/giuseppe/VDR/ALLELESEQ/LC_15/PIPELINE_SANDBOX/d_results_CPo3/interestingHets_NA06986.txt';

if(!$infile){
	print "USAGE: do_alleleseq_remove_interestinghets_in_ENCODEblakclistregions.pl -i=<INFILE>\n";
	print "<INFILE> file interestingHets.txt from AlleleSeq\n";
    exit 1;
}
my ($filename, $directory) = fileparse($infile);
$filename =~ s/(.*)\..*/$1/;
my $temp_bed_file = $directory . $filename . '_temp.bed';
my $temp_out_bed_file = $directory . $filename . 'out_temp.bed'; #here bedtools will put its output
my $out_file = $directory . $filename . '_EBLfiltered.txt';

my $header;
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
open (my $outstream,  q{>}, $temp_bed_file) or die("Unable to open $temp_bed_file : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	if($_ =~ /^chrm/){
		$header = $_;
		next;
	}
	my @fields = split("\t",$_);
	my $start = $fields[1] - 1;
	splice @fields, 1, 0, $start;
	print $outstream join("\t", @fields) . "\n";
}
close $instream;
close $outstream;

#now bed intersect, get the output and rewrite it as it was
system "$BEDTOOLS/bedtools intersect -v -a $temp_bed_file -b $BLACKLIST_REGIONS > $temp_out_bed_file";
open ($instream,  q{<}, $temp_out_bed_file) or die("Unable to open $temp_out_bed_file : $!");
open ($outstream,  q{>}, $out_file) or die("Unable to open $out_file : $!");
print $outstream $header, "\n";
while(<$instream>){
	chomp;
	next if($_ eq '');
	my @fields = split("\t",$_);
	splice @fields, 1, 1;
	print $outstream join("\t", @fields) . "\n";
}
close $instream;
close $outstream;
unlink $temp_bed_file;
unlink $temp_out_bed_file;
