#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#maybe useful to get a bed of asym events to input to gat (although andreas says that gat has no benefit from a hypergeometric test, with snps)

#the alleleseq pipeline outputs a tsv file called "interestingHets.txt"
#this contain the snp position with differential allelic binding
#I want to filter these to obtain all of the following:
#1 snp is within peak and/or at distance d from closest peak
#2 snp is significant (qval <0.05, or whatever the set significance in alleleseq was)
#3 snp is asym


#Ideally I want to pass the ASYM BindingSITE=! snps to SNAP or to Haploreg to pick those in strong LD with them.
#So output a file with chr:pos

#file structured as follows
#chrm    snppos          ref     mat_gtyp        pat_gtyp        c_gtyp  phase   mat_all pat_all cA      cC      cG      cT      winning SymCls  SymPval BindingSite     cnv
#1       797440  T       C       T       Y       HETERO  None    None    0       0       0       8       ?       Asym    0.0078125       0       1.0
#1       1314015 C       Y       C       Y       HETERO  None    None    0       0       0       8       ?       Asym    0.0078125       0       1.0
#1       1371459 A       R       G       R       HETERO  None    None    0       0       6       0       ?       Asym    0.03125 0       1.0
#1       1509156 A       A       G       R       HETERO  None    None    0       0       5       0       ?       Sym     0.0625  0       1.0
#1       1706136 T       Y       T       Y       HETERO  None    None    0       0       0       6       ?       Asym    0.03125 0       1.0

#checks on:
#column 0,1,14,15,16


#The guy from Alleleseq says:

#1) The last column (cnv=1.0) seems characteristic of an 'empty CNV file'
#actually, or at least the dummy value that you inputted manually instead of
#something that is calculated - since usually the normalized read depth
#computation is seldom exactly one. 

#2) Did you feed in a binding site file? Then probably so.

#3) for SymPval, it doesn't mean that they are insignificant; in fact all the SNPs in interestingHEts.txt are supposed to be the "significant" ones. In general, one can perform a simulation, calculate the FDR, for given a certain p-value cutoff. So this was performed in AlleleSeq based on an explicit binomial simulation of the FDR in each AlleleSeq run; there is no single cutoff because of various factors that can influence this, particularly the different sizes of the datasets, number of reads etc. This is actually mentioned in the AlleleSeq paper.  So if you look into your FDR.txt, the very last line says "Target 0.1" (what you instructed before the run), FDR that was cut off (near 0.1) and the corresponding p value that is the threshold. So all the SNPs in the interestingHets.txt file should be below this p threshold.'

my $infile;
GetOptions(
        'i=s'        =>\$infile,
);
if(!$infile){
	print "USAGE: do_alleleseq_get_bed_from_output.pl -i=<INFILE>\n";
	print "<INFILE> file interestingHets.txt from AlleleSeq\n";
    exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
my $counter = 1;
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^chrm/); #header
	
	my ($chr,$snppos,$symCls,$symPval,$bindingSite) = (split /\t/)[0,1,14,15,16];
	next if(!$chr);
	next if(!$snppos);
	next if(!$symCls);
	next if(!$symPval);
	next if(!$bindingSite);
	
	#checks
	next if ($bindingSite ne '1');
	next if ($symCls ne 'Asym');
	
	my $name = 'AlleleSeq_Asym_SNP_' . $counter;
	$chr =~ s/(.+)/chr$1/;
	print $chr . "\t" . ($snppos-1) . "\t" . $snppos . "\t" . $name . "\t" . $symPval . "\n";
	$counter++;  
	
}
close $infile;

