#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#basically just to automate the process of intersecting the VDRucBV (LD clumped) BED SNPs  with the disease LD blocks from GRASP (D' LD)
#inputs
#GRASP PLINK D' LD block beds with disease info
#BED containing LD clumped VDRucBVs and obtained with do_funseq_GAT_VDRhcBVs_to_VDRucBVS.pl

my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";

my $INPUT_VDRBV;
my $INPUT_LD;
GetOptions(
        'i_ld=s'          =>\$INPUT_LD,
        'i_vdrbv=s'       =>\$INPUT_VDRBV
);
if(!$INPUT_LD){
     print "USAGE: do_funseq_intersect_PLINK_LD_GRASP_with_VDRucBVbed.pl -i_ld=<INPUT_LD_FILE> -i_vdrbv=<INPUT_VDRBV>\n";
     print "<INPUT_LD_FILE> gzipped bed containing GRASP catalog snps mapped to D' ld blocks, b37 format\n";
     print "<INPUT_VDRBV> b37 bed file containing VDRucBVs\n";
     exit 1;
}
if(!$INPUT_VDRBV){
     print "USAGE: do_funseq_intersect_PLINK_LD_GRASP_with_VDRucBVbed.pl -i_ld=<INPUT_LD_FILE> -i_vdrbv=<INPUT_VDRBV>\n";
     print "<INPUT_LD_FILE> gzipped bed containing GRASP catalog snps mapped to D' ld blocks, b37 format\n";
     print "<INPUT_VDRBV> b37 bed file containing VDRucBVs\n";
     exit 1;
}
my ($filename,  $directory) = fileparse($INPUT_VDRBV);
my $ethnicity;
if($INPUT_LD =~ /(ceu)/i){
	$ethnicity = $1;
}elsif($INPUT_LD =~ /(yri)/i){
	$ethnicity = $1;
}else{
	print "Error: impossible to determine ethnicity of LD population used from filename: $INPUT_LD.\nAborting..\n";
	exit -1;
}
my $outfile = $directory .  'VDR-ucBVs_INTERSECT_GRASP_LDBLOCKS_'  . $ethnicity . '.tsv';
system "$BEDTOOLS intersect -wo -a $INPUT_LD -b $INPUT_VDRBV  > $outfile";

#bedtools intersect -wo -a LD_PLINK_[CEU|YRI]_GRASP2final.vcf.gz -b interestingHets_consensus_inpeaks_o[1,2,3,4,5,6].vcf.gz > o[1,2,3,4,5,6]_LD_PLINK_[CEU|YRI]_GRASP2final.tsv

system "sed -i '1iCHR	LD_START	LD_STOP	SNPS_IN_BLOCK	BLOCK_SIZE_KB	GRASP_CHR	GRASP_DISEASE_SNP_POS	GRASP_DISEASE_SNP_ID	REF	ALT	QUAL	GRASP_PHENOTYPE	GRASP_PH_CAT_PVAL	BEDTOOLS_INTERSECT	CHR	START	STOP	BEDTOOLS_INTERSECT' $outfile";

print STDERR "Output in $directory\n";
print STDERR "FINISHED.\n";
