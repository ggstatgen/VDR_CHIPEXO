#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#I want a bed version of the SNPS returned by haploreg
#I want to run GAT on these to see how likely they are to intersects features of interest (like GRASP snps, GWAS snps etc)

#input format
#35 query snps
#chr     pos     r2      D'      is_query_snp    rsID    ref     alt     AFR     AMR     ASN     EUR     GERP_cons       SiPhy_cons      Promoter_ENCODE Enhancer_ENCODE #Promoter_Roadmap        Enhancer_Roadmap        DNAse   Proteins
#1       45175925        0.88    0.95    0       rs10890313      G       A       0.43    0.27    0.33    0.15    0       0       .       .       .       .       .       .       .       GR_known3;Hoxa5_1       ENSG00000198520.6       C1orf
#1       45185261        0.92    0.97    0       rs12746898      A       G       0.43    0.27    0.33    0.15    0       0       .       GM12878,6_Weak_Enhancer;HepG2,6_Weak_Enhancer;K562,7_Weak_Enhancer      .       ADI.MSC,11_EnhWk1;CD4  

#coord is hg19 so good

#get, for example, 0,1,4
#and produce a bed version

my $infile;
GetOptions(
        'i=s'        =>\$infile,
);
if(!$infile){
	print "USAGE: do_SNAP_get_bed_from_output.pl -i=<INFILE>\n";
    print "<INFILE> file .txt from HAPLOREG\n";
    exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	next if($_ =~ /query snps/);
	next if($_ =~ /^chr/);
	next if($_ eq '');
	next if($_ =~ /^\s/);
	my ($chr, $pos, $is_query_snp, $rsid) = (split /\t/)[0,1,4,5];
	#build bed
	$chr = 'chr' . $chr;
	print $chr . "\t" . ($pos-1) . "\t" . $pos . "\t" . $rsid . "\t" . $is_query_snp . "\n";
}
close $instream;