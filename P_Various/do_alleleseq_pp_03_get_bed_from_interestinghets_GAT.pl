#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#get a bed of asym events to input to gat (although andreas says that gat has no benefit from a hypergeometric test, with snps)
#give the choice of selecting 1)all peaks 2) in peaks 3) not in peaks

#input is as follows:
#chrm    snppos  ref     BindingSite     overlap sampleIDS
#1       108357419       C       0       2       NA11829,NA19213
#1       11968297        T       1       2       NA19190,NA19213
#1       150459715       T       1       2       NA19189,NA19248
#1       157151799       G       1       4       NA10847,NA12872,NA19247,NA19248
#1       158894077       G       0       2       NA10847,NA11829
#1       167632618       G       1       2       NA19189,NA19213


#bed 
#chr start end name score

my $infile;

GetOptions(
        'i=s'        =>\$infile,
);
if(!$infile){
	print "USAGE: do_alleleseq_pp_03_get_bed_from_interestinghets_GAT.pl -i=<INFILE>\n";
	print "<INFILE> file interestingHets..overlap.txt from do_alleleseq_pp..pl\n";
    exit 1;
}
my ($filename, $directory) = fileparse($infile);
$filename =~ s/(.*)\..*/$1/;
my $outfile_all = $directory . $filename . 'all.bed';
my $outfile_inpeaks  = $directory . $filename . 'inpeaks.bed';
my $outfile_outpeaks = $directory . $filename . 'outpeaks.bed';;

open (my $instream,  q{<}, $infile)                    or die("Unable to open $infile : $!");
open (my $outstream_all,      q{>}, $outfile_all)      or die("Unable to open $outfile_all : $!");
open (my $outstream_inpeaks,  q{>}, $outfile_inpeaks)  or die("Unable to open $outfile_inpeaks : $!");
open (my $outstream_outpeaks, q{>}, $outfile_outpeaks) or die("Unable to open $outfile_outpeaks : $!");
my $counter = 1;
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^chrm/); #header
	
	my ($chr, $snppos, $ref, $bds, $overlap, $sampleids) = (split /\t/)[0,1,2,3,4,5];
	next if(!$chr);
	next if(!$snppos);
	
	my $snp_start = $snppos - 1;
	my $bed_line = 'chr'. $chr . "\t" . $snp_start . "\t" . $snppos . "\t" . $sampleids . "\t" . $overlap;
	
	if($bds =~ /1/){
		print $outstream_all $bed_line, "\n";
		print $outstream_inpeaks $bed_line, "\n";
	}elsif($bds =~ /0/){
		print $outstream_all $bed_line, "\n";
		print $outstream_outpeaks $bed_line, "\n";
	}else{
		print STDERR "Binding site field not recognised: $bds. Skipping..\n";
		next;
	}
}
close $infile;
close $outstream_all;
close $outstream_inpeaks;
close $outstream_outpeaks;