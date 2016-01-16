#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#convert the list of causal autoimmune PICS SNPs found here
#http://www.broadinstitute.org/pubs/finemapping/?q=data-portal
#to bed format. This will then be used as GAT annotation for an overlap analysis

#A simple 4 fields bed should be sufficient: chr start stop disease_name pics_score

#The list has been first converted from xls to tsv
my $INPUT_PICS = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/BROAD_PICS/BROAD_candidate_causal_snps_39immune_nonimmune_diseases_plus_enh_annot_masterfile9.csv";

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



open (my $instream,  q{<}, $INPUT_PICS) or die("Unable to open $INPUT_PICS : $!");
while(<$instream>){
	chomp;
	next if ($_ =~ /^Disease/);
	my @fields = split("\t", $_);
	
	#just recombine the data following the bed format
	#chrom start end name score strand [other]
	my $chr = $fields[3];
	$chr =~ s/chr(.+)/$1/;
	my $stop =  $fields[4];
	my $start = $stop - 1;
	#my $name = $fields[2] . '-' . $fields[0]; #with this GAT wouldn't work
	my $name = $fields[0]; 
	my $score = $fields[5];

	#print $chr . "\t" . $start . "\t" . $stop . "\t" . $name . "\t" . $score . "\n";   
	print $chr . "\t" . $start . "\t" . $stop . "\t" . $name . "\n";                       
}
close $instream;
