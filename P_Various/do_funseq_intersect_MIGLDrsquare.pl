#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

#basically just to automate the process of intersecting the VDR-BV SNPs annotated using funseq with the disease LD blocks for GWAS and GRASP obtained using the MIGLD r2 pipeline
#inputs
#LD block beds with disease info
#VCFs for ASB interesting het snps coming from the allelespecific pipeline

#LD DISEASE INTERVALS
#my $INPUT_GRASP_CEU = "/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/LD_MIGLD_RSQUARE_CEU_GRASP.bed.gz";
my $INPUT_GRASP_CEU = "/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/LD_MIGLD_CEU_r2_0.8_MAF_0.01_frac_0.9_GRASP.bed.gz";
my $input_recur     = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/GWAS_INTERSECTIONS/Output_recur_b37.vcf";

my ($filename,  $directory) = fileparse($input_recur);
my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin/bedtools";

#GRASP data header:
my $header_grasp = "CHR	LD_START	LD_STOP	SNPS_IN_BLOCK	BLOCK_SIZE_KB	GRASP_CHR	GRASP_DISEASE_SNP_POS	GRASP_DISEASE_SNP_ID	REF	ALT	QUAL	GRASP_PHENOTYPE	GRASP_PH_CAT_PVAL	BEDTOOLS_INTERSECT	ASB_CHR	ASB_SNP_POS	ASB_SNP_ID	ASB_REF	ASB_ALT	ASB_QUAL	ASB_FILTER	ASB_1KG_INFO	ASB_1KG_FORMAT	NA06986	NA06989	NA10847	NA11829	NA11919	NA12383	NA12489	NA12872	NA19189	NA19190	NA19213	NA19235	NA19236	NA19247	NA19248	NA07029	NA10831	NA11832	NA19191	NA19249";

#GWAScatalog header:
#my $header_gwas = "CHR	LD_START	LD_STOP	SNPS_IN_BLOCK	BLOCK_SIZE_KB	GWAS_CHR	GWAS_DISEASE_SNP_POS	GWAS_DISEASE_SNP_ID	REF	ALT	QUAL	GWAS_PHENOTYPE	GWAS_PH_CAT_PVAL	BEDTOOLS_INTERSECT	ASB_CHR	ASB_SNP_POS	ASB_SNP_ID	ASB_REF	ASB_ALT	ASB_QUAL	ASB_FILTER	ASB_1KG_INFO	ASB_1KG_FORMAT	NA06986	NA06989	NA10847	NA11829	NA11919	NA12383	NA12489	NA12872	NA19189	NA19190	NA19213	NA19235	NA19236	NA19247	NA19248	NA07029	NA10831	NA11832	NA19191	NA19249";

print STDERR "do_funseq_intersect_MIGLDrsquare:\n";
print STDERR "Using the following files:\n";
print STDERR "GRASP LD blocks CEU: $INPUT_GRASP_CEU\n";
print STDERR "Processing overlap files under $directory\n";
print STDERR "Saving all output under $directory\n";
print STDERR "edit if required\n";

my $out_recur_GRASP_CEU = $directory . "LD_MIGLDrsquare_CEU_GRASP_recur.tsv";
system "$BEDTOOLS intersect -wo -a $INPUT_GRASP_CEU -b $input_recur  > $out_recur_GRASP_CEU";


system "sed -i '1iCHR	LD_START	LD_STOP	SNPS_IN_BLOCK	BLOCK_SIZE_KB	GRASP_CHR	GRASP_DISEASE_SNP_POS	GRASP_DISEASE_SNP_ID	REF	ALT	QUAL	GRASP_PHENOTYPE	GRASP_PH_CAT_PVAL	BEDTOOLS_INTERSECT	ASB_CHR	ASB_SNP_POS	ASB_SNP_ID	ASB_REF	ASB_ALT	ASB_QUAL	ASB_FILTER	ASB_1KG_INFO	ASB_1KG_FORMAT	NA06986	NA06989	NA10847	NA11829	NA11919	NA12383	NA12489	NA12872	NA19189	NA19190	NA19213	NA19235	NA19236	NA19247	NA19248	NA07029	NA10831	NA11832	NA19191	NA19249' $out_recur_GRASP_CEU";

print STDERR "FINISHED.\n";
