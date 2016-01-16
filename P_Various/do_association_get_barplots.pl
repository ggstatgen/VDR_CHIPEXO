#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#I want to get bar plots showing the distribution of normalised read counts from those that SNPTEST tell me are the best candidates

#input 1: TMM normalised matrix of reads
#input 2: snptest output?

#output: some tsv suitable for R/ggplot2
my $infile_rcmatrix;
my $infile_genotypes;
my $infile_rsid;
my $infile_peakid;

my $GENOTYPE_FILE = '/net/isi-scratch/giuseppe/VDR/VARIANTS/IMPUTATION/d_IMPUTE_out_chksize_4mb/IMPUTE_1kg_to_hapmap_autosomes_hg19.vcf.gz';

GetOptions(
        'rc=s'        =>\$infile_rcmatrix,
        'snp=s'       =>\$infile_rsid,
        'peak=s'      =>\$infile_peakid
);

#$infile_rcmatrix = '/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_diffbind/vdr_o3_matrix_p.txt';
#$infile_rsid = 'rs17160772';
#$infile_peakid = '13903';

if(!$infile_rcmatrix){
     print "USAGE: do_association_get_barplots.pl -rc=<RC_MATRIX> -snp=<SNP_RSID> -peak=<PEAK_ID>\n";
     print "<RC_MATRIX> diffbind TMM processed read depths matrix\n";
     print "<SNP_RSID> snp of interest\n";
     print "<PEAK_ID> peak number of interest\n";
     print "Using imputed hg19 genotypes - if not suitable, modify scripts\n";
     exit 1;
}
if(!$infile_rsid){
     print "USAGE: do_association_get_barplots.pl -rc=<RC_MATRIX> -snp=<SNP_RSID> -peak=<PEAK_ID>\n";
     print "<RC_MATRIX> diffbind TMM processed read depths matrix\n";
     print "<SNP_RSID> snp of interest\n";
     print "<PEAK_ID> peak number of interest\n";
     print "Using imputed hg19 genotypes - if not suitable, modify scripts\n";
     exit 1;
}
if(!$infile_peakid){
     print "USAGE: do_association_get_barplots.pl -rc=<RC_MATRIX> -snp=<SNP_RSID> -peak=<PEAK_ID>\n";
     print "<RC_MATRIX> diffbind TMM processed read depths matrix\n";
     print "<SNP_RSID> snp of interest\n";
     print "<PEAK_ID> peak number of interest\n";
     print "Using imputed hg19 genotypes - if not suitable, modify scripts\n";
     exit 1;
}
my %sample_name_to_rd;
my %sample_name_to_gt; 
###############
#read depth processing
###############
my $rc_header = `head -n 1 $infile_rcmatrix`;
chomp $rc_header;
my @sample_names = split("\t", $rc_header);
shift @sample_names;
my $search_string = '"^' . $infile_peakid . '\t"';
my $rd_line = `grep -m 1 -P  $search_string $infile_rcmatrix`;
chomp $rd_line;
my @counts = split("\t", $rd_line);
shift @counts;
@sample_name_to_rd{@sample_names} = @counts;

##############
#Genotype processing
##############
$search_string = '"^\#CHROM"';
my $gt_header = `zgrep -m 1 -P $search_string $GENOTYPE_FILE`; 
chomp $gt_header;
my @cols = split("\t", $gt_header);
my @c_names = @cols[9..(@cols-1)];

#search snp of interest against genotypes
$search_string = '"\W' . $infile_rsid . '\W"';
my $vcf_line = `zgrep -m 1 -P $search_string $GENOTYPE_FILE`; 
chomp $vcf_line;
#map sample names to genotypes
my @vcf_line = split(/\t/, $vcf_line);
my @genotypes = @vcf_line[9..$#vcf_line]; 
#map all samples in the vcf to their genotype:
#(for some of these I won't have normalised reads)
@sample_name_to_gt{@c_names} = @genotypes;

my $snp_chr = $vcf_line[0]; 
my $snp_loc = $vcf_line[1];
my $snp_ref = $vcf_line[3];
my $snp_alt = $vcf_line[4];

my $gt_hr  = $snp_ref . $snp_ref;
my $gt_het = $snp_ref . $snp_alt;
my $gt_hnr = $snp_alt . $snp_alt;
if($snp_alt =~ /(\w\W)+/){
	print "WARNING: alternate allele string: $snp_alt seem to contain multiple alternate alleles. Aborting..\n";
	exit;
}
my $s_to_g = get_sample_to_genotype_hash($snp_ref, $snp_alt, %sample_name_to_gt);

#order by genotype
my %sample_homref_to_line;
my %sample_het_to_line;
my %sample_homnr_to_line;
foreach my $item (sort keys %{$s_to_g}){
	if($$s_to_g{$item} eq $gt_hr){
		$sample_homref_to_line{$item} = $snp_chr . '-' . $infile_rsid  . "\t" . $item . "\t" . $$s_to_g{$item} . "\t" . $sample_name_to_rd{$item};
	}elsif($$s_to_g{$item} eq $gt_het){
		$sample_het_to_line{$item} = $snp_chr . '-' . $infile_rsid . "\t" . $item . "\t" . $$s_to_g{$item} . "\t" . $sample_name_to_rd{$item};
	}elsif($$s_to_g{$item} eq $gt_hnr){
		$sample_homnr_to_line{$item} = $snp_chr . '-' . $infile_rsid . "\t" . $item . "\t" . $$s_to_g{$item} . "\t" . $sample_name_to_rd{$item};
	}elsif($$s_to_g{$item} eq 'NN'){
		print "Warning: sample $item has genotype NN for this SNP. Skipping sample..\n";
		next;
	}
	else{
		print "Error: sample genotype: $sample_name_to_gt{$item} does not match any allele configuration for this snp. Aborting.\n";
	}
		
}

#print according to ordering
foreach my $item (sort keys %sample_homref_to_line){
	print $sample_homref_to_line{$item}, "\n";
}
foreach my $item (sort keys %sample_het_to_line){
	print $sample_het_to_line{$item}, "\n";
}
foreach my $item (sort keys %sample_homnr_to_line){
	print $sample_homnr_to_line{$item}, "\n";
}




###########
#subs
###########
#the desired output is a hash with key the sample name, and value a string as follows:
#X/X
#?? if no data available
#example input
#T C 1/1:0.8888      0/1:0.8888      1/1:0.8888      0/1:0.8888      0/1:0.8004 
sub get_sample_to_genotype_hash{
	my ($my_vcf_ref, $my_vcf_alt, %my_sample_to_gt) = @_;
	my %sample_to_genotype;

	my %allele_map = (
		'0' => $my_vcf_ref,
		'1' => $my_vcf_alt 
	);	

	foreach my $item (sort keys %my_sample_to_gt){
		my @genotype_fields = split(':', $my_sample_to_gt{$item});
		my $sample_gt = $genotype_fields[0];
		#check that the genotype is actually present: 
		if($sample_gt =~ /(\d{1})[\/|\|](\d{1})/){
			if(!$allele_map{$1}) { print "Error: genotype $sample_gt contains more than one alternate allele: $1. Aborting..\n"; }
			if(!$allele_map{$2}) { print "Error: genotype $sample_gt contains more than one alternate allele: $2. Aborting..\n"; }
			my $output = $allele_map{$1}  .  $allele_map{$2};
			$sample_to_genotype{$item} = $output;
		}else{
			#genotype no present
			$sample_to_genotype{$item} = 'NN';
		}
	}
	return \%sample_to_genotype;
}
