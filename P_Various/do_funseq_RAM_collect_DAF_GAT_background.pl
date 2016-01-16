#!/usr/bin/perl
use strict;
use warnings;
use IO::Zlib;
use List::MoreUtils qw(:all);

#ALSO get GLOBAL MAF

#Chris email, "VDR paper" 04/12/15 18:56
#This accounts for the pre-existing enrichment in similar annotation observed at VDR binding locations and additionally he took account of GC-composition effects. 
#Thinking more about this, it's not entirely clear to me why a minor/derived allele frequency-matched background is the most appropriate here, and Fahr et al.'s (frequency-based) approach to the background is not as sophisticated as the one that Giuseppe now uses. 

#####
#We could post-hoc check that the distribution of MAFs for all 1000G SNPs under ChIP-Exo peaks, and the distribution of MAFs for all VDR-BVs, are drawn from the same distribution (KS test). (The GAT-based enrichment test we employed is not coded to take account of frequency variations, so making a change at this point would put the work back substantially and we could not submit in the next month as we now wish to do.)
#####

#This script will collect
#0 MAF for all non 1000g SNPs in the 20 samples used for alleleseq (13,151,808) 
#1 MAFs for all 1000g in the background I used. Currently, this is: VDR-ASB background (sym/asym,>=5 reads) and VDR-QTL background 
#2 MAFs for all VDR-hcBVs (1,695 positions)
#3 MAFs for all VDR (~43,000 positions)

#UNFORTUNATELY, I have some imputed variables and don't have 1kg frequencies for them.

#Then call R and do a comparison of the distributions

#info from 1000genomes
##INFO=<ID=AC,Number=.,Type=Integer,Description="Alternate Allele Count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total Allele Count">
##INFO=<ID=AF,Number=1,Type=Float,Description="Global Allele Frequency based on AC/AN">
##INFO=<ID=AMR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from AMR based on AC/AN">
##INFO=<ID=ASN_AF,Number=1,Type=Float,Description="Allele Frequency for samples from ASN based on AC/AN">
##INFO=<ID=AFR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from AFR based on AC/AN">
##INFO=<ID=EUR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from EUR based on AC/AN">

#ATTENTION: when the ancestral allele is the alternate allele, the info above is the ANCESTRAL ALLELE COUNT.
#To get the DERIVED ALLELE COUNT, do a 1-ANCESTRAL ALLELE COUNT

#STRATEGY: did a bedtool intersect between
#1 1000g variants tested AND BACKGROUND -> vcf1
#2 1000g variants tested and VDR-BVs -> vcf2
#for details see README in /net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES

my $DIRECTORY = '/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/d_GAT_BACKGROUNDS/';

#This is a simple script to collect MAF from these two files
my $INPUT_VDRHCBV             = $DIRECTORY . 'VDR-hcBV_1kgdata_b37_1kg_data.vcf.gz'; #b37, 1695 intervals
my $INPUT_VDRBV               = $DIRECTORY . 'VDR-BV_1kgdata_b37_1kg_data.vcf.gz'; #b37 43332 intervals
my $INPUT_BACKGROUND_FIG5_A   = $DIRECTORY . 'GAT_BKG_FIG5A_ALLELESEQ_20samples_merged_INTERSECT_cpo3.vcf.gz'; #b37, 20,356 variants in cpo3 peaks
my $INPUT_BACKGROUND_FIG5_BCD = $DIRECTORY . "GAT_BACKGROUND_SIMASYM_FINAL.vcf.gz"; #b37, 114,142 variants (all sym/asym, >= 5 reads pileup)
#my $INPUT_BACKGROUND_SYMASYM = $DIRECTORY . 'GAT_BACKGROUND_SIMASYM_QTL_1kg_noindels.vcf.gz'; #b37, 159,413 variants (all tested asym/sym, all tested BIMBAM/SNPTEST)
#my $INPUT_BACKGROUND_ALL = "/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/GAT_BACKGROUND_ALL.vcf.gz";  #b37, 114,142 variants (all sym/asym, >= 5 reads pileup)

#sample lines
#1       618463  rs114496243     G       A       100     PASS    AA=.;AF=0.2;AFR_AF=0.34;AMR_AF=0.17;ASN_AF=0.24;AVGPOST=0.9066;ERATE=0.0107;EUR_AF=0.09;LDAF=0.2272;RSQ=0.7682;SNPSOURCE=LOWCOV;THETA=0.0285;VT=SNP;AN=30;AC=6
#1       714019  rs114983708     A       G       100     PASS    AA=.;AF=0.12;AFR_AF=0.37;AMR_AF=0.05;ASN_AF=0.04;AVGPOST=0.983;ERATE=0.0013;EUR_AF=0.04;LDAF=0.1171;RSQ=0.9494;SNPSOURCE=LOWCOV;THETA=0.0038;VT=SNP;AN=30;AC=4

#Open each of these and create and put results in R file

#before saving the frequencies, you need to to KNOW if the ancestral allele is the ref or the alt
#if the ancestral allele is the ref, save the frequency as is
#if the ancestral allele is the alt, the frequency you have is for the ancestral. The derived will be (1 - freq)



#############
#DAF, VDR-hcBV
#############
my @variants_vdrhcbv_eur; my @variants_vdrhcbv_afr;my @variants_vdrhcbv_glob;
get_daf_from_vcf_variants($INPUT_VDRHCBV, \@variants_vdrhcbv_eur, \@variants_vdrhcbv_afr, \@variants_vdrhcbv_glob);
#############
#DAF, VDR-BV
#############
my @variants_vdrbv_eur; my @variants_vdrbv_afr;my @variants_vdrbv_glob;
get_daf_from_vcf_variants($INPUT_VDRBV, \@variants_vdrbv_eur, \@variants_vdrbv_afr, \@variants_vdrbv_glob);
#############
#DAF, VDR Background, Fig 5A (~20,000 variants)
#############
my @variants_vdrbkg_fig5a_eur; my @variants_vdrbkg_fig5a_afr;my @variants_vdrbkg_fig5a_glob;
get_daf_from_vcf_variants($INPUT_BACKGROUND_FIG5_A, \@variants_vdrbkg_fig5a_eur, \@variants_vdrbkg_fig5a_afr, \@variants_vdrbkg_fig5a_glob);
#############
#DAF, VDR Background, Fig 5B,C,D (~114,000 variants)
#############
my @variants_vdrbkg_fig5bcd_eur; my @variants_vdrbkg_fig5bcd_afr;my @variants_vdrbkg_fig5bcd_glob;
get_daf_from_vcf_variants($INPUT_BACKGROUND_FIG5_BCD, \@variants_vdrbkg_fig5bcd_eur, \@variants_vdrbkg_fig5bcd_afr, \@variants_vdrbkg_fig5bcd_glob);
############
#DAF, All 1kg Variants used in Alleleseq Test (millions of varianta, includes variants in coding regions and elsewhere which have no ChIP-exo coverage)
############
#my @variants_bkg_all_eur;my @variants_bkg_all_afr;
#get_daf_from_vcf_variants($INPUT_BACKGROUND_ALL, \@variants_bkg_all_eur, \@variants_bkg_all_afr);


#choose right header
#print "ETHNICITY\tDAF_VDRhcBV\tDAF_VDRBV\tDAF_VDRBKG\tDAF_1KGALL\n";
#print "ETHNICITY\tDAF_VDRhcBV\tDAF_VDRBV\tDAF_VDRBKG\n";
#print "ETHNICITY\tDAF_VDRhcBV\tDAF_VDRBV\tDAF_VDRBKG_5A\tDAF_VDRBKG_5BCD\n";
print "ETHNICITY\tDAF_VDRhcBV\tDAF_VDRBKG_5BCD\n"; #ETHNICITY = CEU, YRI, GLOBAL


#create tsv
#excellent solution found here:
#http://www.perlmonks.org/?node_id=446821
my $result_ceu     = each_array(@variants_vdrhcbv_eur,  @variants_vdrbv_eur,  @variants_vdrbkg_fig5a_eur,  @variants_vdrbkg_fig5bcd_eur);
my $result_yri     = each_array(@variants_vdrhcbv_afr,  @variants_vdrbv_afr,  @variants_vdrbkg_fig5a_afr,  @variants_vdrbkg_fig5bcd_afr);
my $result_glob    = each_array(@variants_vdrhcbv_glob, @variants_vdrbv_glob, @variants_vdrbkg_fig5a_glob, @variants_vdrbkg_fig5bcd_glob);
while ( my ($a, $b, $c, $d) = $result_ceu->() ){

        $a = 'NA' unless($a);
        $b = 'NA' unless($b);
        $c = 'NA' unless($c);
        $d = 'NA' unless($d);

    #print "CEU\t$a\t$b\t$c\t$d\n";
    print "CEU\t$a\t$d\n";
}
while ( my ($a, $b, $c, $d) = $result_yri->() ){

        $a = 'NA' unless($a);
        $b = 'NA' unless($b);
        $c = 'NA' unless($c);
        $d = 'NA' unless($d);

    #print "YRI\t$a\t$b\t$c\t$d\n";
    print "YRI\t$a\t$d\n";
}
while ( my ($a, $b, $c, $d) = $result_glob->() ){

        $a = 'NA' unless($a);
        $b = 'NA' unless($b);
        $c = 'NA' unless($c);
        $d = 'NA' unless($d);

    #print "GLOB\t$a\t$b\t$c\t$d\n";
    print "GLOB\t$a\t$d\n";
}

#################################################################


sub get_daf_from_vcf_variants{
	my ($INPUT_FILE, $array_eur_ref, $array_afr_ref, $array_glob_ref) = @_;

	tie *FILE,   'IO::Zlib', $INPUT_FILE, "rb";
	while (<FILE>)	{ 
		chomp;
		next if($_ eq '');
		next if($_ =~ /^#/);
		my $ALTF_AFR; my $ALTF_EUR; my $ALTF_GLOB;
		my $FLAG; # set to one if the ancestral is the alternate	
	
		my @fields = split("\t", $_);
		#skip indels
		next if(length($fields[3]) > 1);
		next if(length($fields[4]) > 1);
		my $ref = $fields[3]; my $alt = $fields[4]; my $info = $fields[7];
	
		#get ancestral allele info=====================
		#if($info[0] =~ /AA=(.*);/){
		if($fields[7] =~ /AA=(\w+)/){			
			my $anc = uc($1);
			if($anc eq $alt){
				$FLAG = 1;		
			}elsif($anc eq $ref){			
			}elsif( ($anc eq '') or ($anc eq 'N') or ($anc eq '.') or ($anc eq '-') ){
				next;
			}else{ 
				next;			
			}
		}else{
			next;		
		}
		#get allele frequencies========================
		my @info = split(";", $fields[7]);
		foreach my $item (@info){
			if($item =~ /^AFR_AF=([0-9]+\.[0-9]+)/){
				$ALTF_AFR = $1;				
			}	
			if($item =~ /^EUR_AF=([0-9]+\.[0-9]+)/){
				$ALTF_EUR = $1;	
			}
			if($item =~ /^AF=([0-9]+\.[0-9]+)/){
				$ALTF_GLOB = $1;
			}
		}	
		
		if($ALTF_EUR){
			if($FLAG){
				push(@$array_eur_ref, (1 - $ALTF_EUR));			
			}else{
				push(@$array_eur_ref, $ALTF_EUR);
			}
		}
		if($ALTF_AFR){
			if($FLAG){
				push(@$array_afr_ref, (1 - $ALTF_AFR));	
			}else{
				push(@$array_afr_ref, $ALTF_AFR);
			}
		}
		if($ALTF_GLOB){
			if($FLAG){
				push(@$array_glob_ref, (1 - $ALTF_GLOB));	
			}else{
				push(@$array_glob_ref, $ALTF_GLOB);
			}
		}		
	}
	close FILE;	
}