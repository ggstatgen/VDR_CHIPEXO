#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#basically just to automate the process of intersecting the VDR-BV SNPs annotated using funseq with the disease LD blocks for GWAS and GRASP
#inputs
#LD block beds with disease info
#VCFs for ASB interesting het snps coming from the allelespecific pipeline

my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";

#GRASP data header:
my $header_grasp = "CHR	LD_START	LD_STOP	SNPS_IN_BLOCK	BLOCK_SIZE_KB	GRASP_CHR	GRASP_DISEASE_SNP_POS	GRASP_DISEASE_SNP_ID	REF	ALT	QUAL	GRASP_PHENOTYPE	GRASP_PH_CAT_PVAL	BEDTOOLS_INTERSECT	ASB_CHR	ASB_SNP_POS	ASB_SNP_ID	ASB_REF	ASB_ALT	ASB_QUAL	ASB_FILTER	ASB_1KG_INFO	ASB_1KG_FORMAT	NA06986	NA06989	NA10847	NA11829	NA11919	NA12383	NA12489	NA12872	NA19189	NA19190	NA19213	NA19235	NA19236	NA19247	NA19248	NA07029	NA10831	NA11832	NA19191	NA19249";


my $INPUT_VDRBV;
my $INPUT_LD;
GetOptions(
        'i_ld=s'          =>\$INPUT_LD,
        'i_vdrbv=s'       =>\$INPUT_VDRBV
);
if(!$INPUT_LD){
     print "USAGE: do_funseq_intersect_PLINK_LD_GRASP.pl -i_ld=<INPUT_LD_FILE> -i_vdrbv=<INPUT_VDRBV>\n";
     print "<INPUT_LD_FILE> gzipped bed containing GRASP catalog snps mapped to D' ld blocks, b37 format\n";
     print "<INPUT_VDRBV> vcf file containing Funseq-annotated VDR-BVs b37 format\n";
     exit 1;
}
if(!$INPUT_VDRBV){
     print "USAGE: do_funseq_intersect_PLINK_LD_GRASP.pl -i_ld=<INPUT_LD_FILE> -i_vdrbv=<INPUT_VDRBV>\n";
     print "<INPUT_LD_FILE> gzipped bed containing GRASP catalog snps mapped to D' ld blocks, b37 format\n";
     print "<INPUT_VDRBV> vcf file containing Funseq-annotated VDR-BVs\n";
     exit 1;
}
#LD DISEASE INTERVALS
#my $INPUT_GRASP_CEU = "/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/LD_PLINK_CEU_GRASP2final.vcf.gz";
#my $INPUT_GRASP_YRI = "/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/LD_PLINK_YRI_GRASP2final.vcf.gz";
#my $input_all       = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/GWAS_INTERSECTIONS/Output_b37.vcf";
#my $input_vdr       = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/GWAS_INTERSECTIONS/Output_invdrpeaks_b37.vcf";
#my $input_recur     = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/GWAS_INTERSECTIONS/Output_recur_b37.vcf";
#my $input_vdr_recur = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/GWAS_INTERSECTIONS/Output_recur_invdrpeaks_b37.vcf";
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
my $out_recur        = $directory .  'LD_PLINK_GRASP_'           . $ethnicity . '_recur.tsv';
#my $out_all_GRASP_CEU = $directory . "LD_PLINK_CEU_GRASP_all.tsv";
#my $out_vdr_GRASP_CEU = $directory . "LD_PLINK_CEU_GRASP_vdr.tsv";
system "$BEDTOOLS intersect -wo -a $INPUT_LD -b $INPUT_VDRBV  > $out_recur";

#bedtools intersect -wo -a LD_PLINK_[CEU|YRI]_GRASP2final.vcf.gz -b interestingHets_consensus_inpeaks_o[1,2,3,4,5,6].vcf.gz > o[1,2,3,4,5,6]_LD_PLINK_[CEU|YRI]_GRASP2final.tsv

system "sed -i '1iCHR	LD_START	LD_STOP	SNPS_IN_BLOCK	BLOCK_SIZE_KB	GRASP_CHR	GRASP_DISEASE_SNP_POS	GRASP_DISEASE_SNP_ID	REF	ALT	QUAL	GRASP_PHENOTYPE	GRASP_PH_CAT_PVAL	BEDTOOLS_INTERSECT	ASB_CHR	ASB_SNP_POS	ASB_SNP_ID	ASB_REF	ASB_ALT	ASB_QUAL	ASB_FILTER	ASB_1KG_INFO	ASB_1KG_FORMAT	NA06986	NA06989	NA10847	NA11829	NA11919	NA12383	NA12489	NA12872	NA19189	NA19190	NA19213	NA19235	NA19236	NA19247	NA19248	NA07029	NA10831	NA11832	NA19191	NA19249' $out_recur";

print STDERR "Output in $directory\n";
print STDERR "FINISHED.\n";
