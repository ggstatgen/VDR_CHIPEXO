#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#Postprocessing of alleleseq output.
#Given a set of alleleseq asym SNP, you want to do a motif analysis of allele specific VDR sites. This is used to determine the effects of the SNPs on the VDR:RXR motif
#1 For all the SNPS, generate 2 sequences containing either SNP and a window centered on the SNP (using bedtools)
#2 analyse these sequences for the presence of the VDR:RXR motif  
#3 then, cross-correlate the difference in the motif score between the two alleles with the degree of allele-specific bias for each snp

#This script outputs the fasta files

#idea taken from supplemental figure 10 in 
#http://www.sciencemag.org/content/suppl/2010/03/18/science.1184655.DC1/McDaniell.SOM.pdf

#INPUTS
#paternal fasta
#maternal fasta
#interesting hets from alleleseq

#TOOLS
#FIMO local?
#bedtools
my $BEDTOOLS = '/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin';
#my $CHROM_SIZE_PATH = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes';
my $CHROM_SIZE_PATH = '/net/isi-scratch/giuseppe/indexes/chrominfo/g1k_v37_chrom.sizes';


#USE DO_FASTA_RENAME_HEADERS.pl on the following two files if you have just cat'ed them from the vcf2diploid tool output
#NA19249
my $PATERNAL_GENOME_PATH = '/net/isi-scratch/giuseppe/VDR/VARIANTS/IMPUTATION/d_IMPUTE_out_chksize_4mb/d_ALLELESEQ/d_VCF2DIPLOID_g1k_v37/YRI_Y120/NA19249_paternal.fa';
my $MATERNAL_GENOME_PATH = '/net/isi-scratch/giuseppe/VDR/VARIANTS/IMPUTATION/d_IMPUTE_out_chksize_4mb/d_ALLELESEQ/d_VCF2DIPLOID_g1k_v37/YRI_Y120/NA19249_maternal.fa';

#NA19191
#my $pat_basename = "_NA19191_paternal.fa";
#my $mat_basename = "_NA19191_maternal.fa";
#my $PERSONAL_GENOME_PATH="/net/isi-scratch/giuseppe/VDR/VARIANTS/omni25_hg19_b37/d_TRIO_DATA_b37/d_VCF2DIPLOID_g1k_v37_1_XYM/YRI_Y111";

#the alleleseq pipeline outputs a tsv file called "interestingHets.txt"
#file structured as follows
#chrm    snppos          ref     mat_gtyp        pat_gtyp        c_gtyp  phase   mat_all pat_all cA      cC      cG      cT      winning SymCls  SymPval BindingSite     cnv
#1       797440  T       C       T       Y       HETERO  None    None    0       0       0       8       ?       Asym    0.0078125       0       1.0
#1       1314015 C       Y       C       Y       HETERO  None    None    0       0       0       8       ?       Asym    0.0078125       0       1.0
#1       1371459 A       R       G       R       HETERO  None    None    0       0       6       0       ?       Asym    0.03125 0       1.0
#1       1509156 A       A       G       R       HETERO  None    None    0       0       5       0       ?       Sym     0.0625  0       1.0
#1       1706136 T       Y       T       Y       HETERO  None    None    0       0       0       6       ?       Asym    0.03125 0       1.0

#will need
#column 0,1,14,15

my $in_asymsnps;
my $SLOP;
GetOptions(
        'snps=s'     =>\$in_asymsnps,
        'slop=i'     =>\$SLOP
);
if(!$in_asymsnps){
	print "USAGE: do_alleleseq_snp_FIMO.motif.pl -snps=<ASYM_SNPS> -slop=<BASEPAIRS>\n";
	print "<ASYM_SNPS> file interestingHets.txt from AlleleSeq\n";
    print "OPTIONAL - <BASEPAIRS> number of bp requested upstream (downstream) of each SNP (default: 20bp+/-)\n";
    exit 1;
}
$SLOP = 20 if(!$SLOP);
my ($filename, $directory) = fileparse($in_asymsnps);
$filename =~ s/(.*)\..*/$1/;
my $out_bed = $directory . "\/" .  $filename . '_temp.bed';
my $out_fasta_pat = $directory . "\/" .  $filename . '_pat.fa';
my $out_fasta_mat = $directory . "\/" .  $filename . '_mat.fa';

#get all good snps from the alleleseq output, create a bed, slop it, get 2 fasta sequences
#my %async_bed;
open (my $instream,  q{<}, $in_asymsnps) or die("Unable to open $in_asymsnps : $!");
open (my $tempstream,  q{>}, $out_bed) or die("Unable to open $out_bed : $!");
while(<$instream>){
	chomp;
	my ($chr,$snppos,$symCls,$symPval,$bindingSite) = (split /\t/)[0,1,14,15,16];
	next if(!$chr);
	next if(!$snppos);
	next if(!$symCls);
	next if(!$symPval);
	next if ($symCls ne 'Asym');
	next unless ($bindingSite =~ /1/);
	
	my $snp_interval =  $chr . "\t" . ($snppos-1) . "\t" . $snppos;
	print $tempstream $snp_interval, "\n";
	#$async_bed{$snp_interval} = 1;
	
	#using process substitution, for future reference
	#http://tldp.org/LDP/abs/html/process-sub.html
	#system "echo -e $snp_interval | $BEDTOOLS/bedtools slop -g $CHROM_SIZE_PATH -b $SLOP -i stdin | $BEDTOOLS/bedtools getfasta -fi $paternal_fasta_temp -fo ";
}
close $instream;
close $tempstream;

#open (my $tempstream,  q{>}, $out_bed) or die("Unable to open $out_bed : $!");
#foreach my $item (sort keys %async_bed){
	#print $tempstream $item, "\n";
#}
#close $tempstream;

#get fasta from paternal and maternal genomes
system "$BEDTOOLS/bedtools slop -g $CHROM_SIZE_PATH -b $SLOP -i $out_bed | $BEDTOOLS/bedtools getfasta -fi $PATERNAL_GENOME_PATH -bed stdin -fo $out_fasta_pat";
system "$BEDTOOLS/bedtools slop -g $CHROM_SIZE_PATH -b $SLOP -i $out_bed | $BEDTOOLS/bedtools getfasta -fi $MATERNAL_GENOME_PATH -bed stdin -fo $out_fasta_mat";

