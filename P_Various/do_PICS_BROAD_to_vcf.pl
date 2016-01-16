#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#convert the list of causal autoimmune PICS SNPs found here
#http://www.broadinstitute.org/pubs/finemapping/?q=data-portal
#to vcf format for intersection with your ASB + QTL snps.

my $INPUT_PICS = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/BROAD_PICS/BROAD_candidate_causal_snps_39immune_nonimmune_diseases_plus_enh_annot_masterfile9.csv";

#in a minimal vcf to be used with bedtools (I want to intersect it with my peaks)
#The MINIMUM vcf understood by bedtools contains the following fields
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

#The list has been first converted from xls to tsv
#fields:

#Disease
#IndexSNP_riskAllele
#SNP
#chr
#pos
#PICS_probability
#Annotation
#nearestGene
#eQTL
#eQTLdirection
#topEnhancer
#inferior_temporal_lobe
#angular_gyrus
#mid_frontal_lobe
#cingulate_gyrus
#substantia_nigra
#anterior_caudate
#hippocampus_middle
#CD25-_CD45RA+_naive
#CD25-_CD45RO+_mem
#CD25+_CD127-_Treg
#CD25-_IL17-_Th_stim_MACS
#CD25-_IL17+_Th17_stim
#Th0
#Th1
#Th2
#CD45RA_CD8
#CD45RO_CD8
#adult_CD14
#adult_CD20
#GM12878
#B_Cell_Centroblast
#Mobilized_CD34
#K562
#colonic_mucosa
#Duodenum_mucosa
#AdiposeNuclei
#HepG2
#Liver1
#Pancreatic_Islets
#kidney
#HSMM-myotube
#NH-Osteoblast
#chondrogenic_dif_cells

my $INFO_DISEASE_def = 'DISEASE='; #field 0
my $INFO_RISK_ALLELE_def = 'RISK_ALLELE='; #field 1 - modify
my $INFO_ANNOTATION_def = 'ANNOTATION='; #field 6
my $INFO_NEAREST_GENE_def = 'NEAREST_GENE='; #field 7
my $INFO_EQTL_def = 'EQTL='; #field 8
my $INFO_EQTL_DIR_def = 'EQTL_DIR='; #field 9
my $INFO_TOP_ENH_def = 'TOP_ENH='; #field 10
my $INFO_GM12878_def = 'GM12878'; #field 30; bool;

#print header
print "##fileformat=VCFv4.1\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";



open (my $instream,  q{<}, $INPUT_PICS) or die("Unable to open $INPUT_PICS : $!");
while(<$instream>){
	chomp;
	next if ($_ =~ /^Disease/);
	my @fields = split("\t", $_);
	
	#just recombine the data following the bed format
	#chrom start end name score strand [other]
	my $chr = $fields[3];
	
	#$chr =~ s/chr(.+)/$1/;
	
	my $pos =  $fields[4];
	my $rsID = $fields[2];
	
	#INFO bits of data
	#get risk allele
	my $risk_allele;
	if($fields[1] =~ /^rs(\d+)-([ACGT])/){
		$risk_allele = $2;
	}else{
		$risk_allele = '-';
	}
	
	my $INFO_RISK_ALLELE = $INFO_RISK_ALLELE_def . $risk_allele; 
	my $INFO_ANNOTATION = $INFO_ANNOTATION_def . $fields[6]; #if none should it not be there?
	my $INFO_NEAREST_GENE = $INFO_NEAREST_GENE_def . $fields[7]; #if none, should it not be there
	my $INFO_EQTL = $INFO_EQTL_def . $fields[8]; #if none, should it not be there
	my $INFO_EQTL_DIR = $INFO_EQTL_DIR_def . $fields[9]; #if none, should it not be there
	my $INFO_TOP_ENH = $INFO_TOP_ENH_def . $fields[10]; #if none, should it not be there?
	my $bool_gm12878 = $fields[30];

	my @INFO_ARRAY;
	push @INFO_ARRAY, $INFO_RISK_ALLELE,$INFO_ANNOTATION,$INFO_NEAREST_GENE,$INFO_EQTL,$INFO_EQTL_DIR,$INFO_TOP_ENH;
	push @INFO_ARRAY, $INFO_GM12878_def if($bool_gm12878 eq '1');

	my $INFO = join(";", @INFO_ARRAY); 

	#build data line
	#currently REF and ALT are set to '-'
	#QUAL will contain the PICS score (field 5)
	#FILTER will contain the DISEASE (field 0)
	my $REF = '-';
	my $ALT = '-';
	my $QUAL = $fields[5];
	my $FILTER = $fields[0];
	
	my $data_line = $chr  . "\t" .
					$pos  . "\t" . 
					$rsID . "\t" . 
					$REF  . "\t" .  
	                $ALT  . "\t" .
	                $QUAL . "\t" . 
	                $FILTER   . "\t" .
	                $INFO;
	print $data_line, "\n";	                          
}
close $instream;
