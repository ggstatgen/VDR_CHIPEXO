#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#bimbam for bayes factor based SNP-phenotype association
#https://www.bcm.edu/research/centers/childrens-nutrition-research-center/mcmcmc/index.cfm?pmid=18981

#This script produces the required input genotype file
#eg
#5 individuals, 4 snps
#5
#4
#IND, id1, id2, id3, id4, id5
#rs1, AT, TT, ??, AT, AA
#rs2, GG, CC, GG, CC, GG
#...

#and the required SNP location file
#rs1, 1200, 1
#where the third is the CHROMOSOME NUMBER

#so you will need
#for the snp location file:
#columns 0,1,2
#for the genotype file
#columns 1+genotypes

#input genotypes
#these have been obtained as follows
#bedtools intersect 
#-a <GENOTYPES.vcf>  
#-b /net/isi-scratch/giuseppe/VDR/POOLED_4/08_BAM_STAMPY_g1k_v37_MIN21/d_gps_smooth3_mrc20/d_bed_RBL/d_diffbind/vdr_o3_CONSENSUS.bed
#so it represents all the Omni 2.5 snps over the available 23 samples (respectively, the hapmap snps over the available 30 samples), intersecting with a consensus peakset obtained with diffbind with min_overlap = 3 samples/peak

#INPUT 1: omni 2.5 vcf
#my $input_vcf = "/net/isi-scratch/giuseppe/VDR/VARIANTS/omni25_hg19_b37/Omni25_genotypes_2141_samples.hg19_OMNI_25_GEN_PASS_INTERSECT_vdr_o3_CONSENSUS.vcf";
#INPUT 2 : hapmap r27 vcf
#my $input_vcf = "/net/isi-scratch/giuseppe/VDR/VARIANTS/hapmap_liftover_hg19_b37/genotypes_CEUYRI_30_r27_nr.hg19_fwd_INTERSECT_vdr_o3_CONSENSUS.vcf";
#INPUT 3 : 1000g
#my $input_vcf = "/net/isi-scratch/giuseppe/VDR/VARIANTS/1000g_hg19_b37/ALL.phase1_release_v3.20101123.snps_indels_svs.genotypes_FIFTEEN_1K_GEN_hg19_INTERSECT_vdr_o3_CONSENSUS.vcf";


#NOTE I'm removing larger variation, keeping only SNP (does bimbam work with larger indels?)

#based on the results of the population analysis, there are some samples which are suspicious (low correlation with other samples, low number of mapping reads, etc)
#here I create a set of samples which are bad, and won't appear in the phenotype counts - so they should not appear in the genotype counts
#you will have to change this set based on the samples you want to discard
my %DISCARDED_SAMPLESET = (
	'10846' => 1,
	'11919' => 1,
	'06997' => 1,
	'12752' => 1
);

my $input_vcf;
my $outdir;
GetOptions(
		'vcf=s'            =>\$input_vcf,
        'outdir=s'         =>\$outdir
);
if(!$input_vcf){
     print "USAGE: do_DistiLD_to_bed.pl -vcf=<INPUT_VCF> -outdir=<OUTDIR>\n";
     print "<INPUT_VCF> input vcf file\n";
     print "<OUTDIR> output directory\n";
     exit 1;
}
if(!$outdir){
     print "USAGE: do_DistiLD_to_bed.pl -vcf=<INPUT_VCF> -outdir=<OUTDIR>\n";
     print "<INPUT_VCF> input vcf file\n";
     print "<OUTDIR> output directory\n";
     exit 1;
}
my($basename, $directory) = fileparse($input_vcf);
$basename =~ s/(.*)\..*/$1/;
my $outfile_gt      = $outdir . "\/" . 'bimbam_'    . $basename . '_gt.txt';
my $outfile_snploc  = $outdir . "\/" . 'bimbam_'    . $basename . '_snploc.txt';

my $snp_counter;
my @cols;my @c_names; my $cell_number; my $sample_names;
my %sample_to_gt;
my %snp_collection;
#open outputs as well
open (my $instream,     q{<}, $input_vcf) or die("Unable to open $input_vcf : $!");
while(<$instream>){
	chomp;
	next if ($_ =~ /^\#\#/); #metadata

	if($_ =~ /^\#CHROM/){ #header
			@cols = split("\t", $_);
		    @c_names = @cols[9..(@cols-1)];
		    next;
	}
	
	my @snp_line = split(/\t/, $_);
	#collect coordinates
	my $snp_chr = $snp_line[0];
	next if(!$snp_chr);
	my $snp_pos = $snp_line[1];
	next if(!$snp_pos);
	my $snp_id  = $snp_line[2];
	next if(!$snp_id);
	next if ($snp_id eq '.');
	my $snp_ref  = $snp_line[3];
	next if(!$snp_ref);
	next if($snp_ref =~ /^[ACGT]{2,}/);
	my $snp_alt  = $snp_line[4];
	next if(!$snp_alt);
	next if($snp_alt =~ /^[ACGT]{2,}/);
	
	my $snp_info = $snp_line[7];
	my $snp_ancestral = get_ancestral_allele($snp_info); # NOT USED see get_sample_to_genotype_hash
	
	#collect genotypes
	my @genotypes = @snp_line[9..(@cols-1)]; 
	@sample_to_gt{@c_names} = @genotypes;
	
	#the next function assumes there is ONLY ONE ALTERNATE ALLELE (not a list)
	#if this is not the case, next
	if($snp_alt =~ /(\w\W)+/){
		print "ERROR: alternate allele string: $snp_alt seem to contain multiple alternate alleles. Skipping..\n";
		next;
	}
	my $s_to_g = get_sample_to_genotype_hash($snp_ref, $snp_alt, %sample_to_gt);
	$snp_collection{$snp_id}{'chr'} =  $snp_chr;
	$snp_collection{$snp_id}{'pos'} =  $snp_pos;
	foreach my $sample (sort keys % {$s_to_g} ) {
		$snp_collection{$snp_id}{'gt'}{$sample} = $$s_to_g{$sample};
	}
	#$snp_collection{$snp_id}{'gts'} =  $gts;
	#print or save the line in a meta hash
	$snp_counter++;
}
close $instream;

####
#genotype file
####
open (my $out_gt_stream,  q{>}, $outfile_gt) or die("Unable to open $outfile_gt : $!");
#header:  number of individuals, number of snps, individuals names
$cell_number = @cols - 9;
$sample_names = join(", ", @c_names);
print $out_gt_stream $cell_number, "\n";
print $out_gt_stream $snp_counter, "\n";
#take it from the hash?
print $out_gt_stream 'IND, ' . $sample_names . "\n";
#print data
foreach my $this_snp_id ( sort keys %snp_collection ){
	print $out_gt_stream $this_snp_id;
	foreach my $this_sample ( sort keys %{ $snp_collection{$this_snp_id}{'gt'} } ){
		print $out_gt_stream ', ' . $snp_collection{$this_snp_id}{'gt'}{$this_sample};
	}
	print $out_gt_stream "\n";
}
close $out_gt_stream;

####
#snp location file
####
open (my $out_snploc_stream,  q{>}, $outfile_snploc) or die("Unable to open $outfile_snploc : $!");
foreach my $this_snp_id ( sort keys %snp_collection ){
	#turn chromosome into number
	my $chr = $snp_collection{$this_snp_id}{'chr'};
	if($chr =~ /^chr(\d+)/){
		$chr = $1;
	}elsif($chr =~ /^chrX/i){
		$chr = '23';
	}elsif($chr =~ /^chrY/i){
		$chr = '24';
	}elsif($chr =~ /^chrM/i){
		$chr = '25';
	}else{
		print "Unrecognised chromosome string: $chr. Aborting.\n";
		exit;
	}
	
	print $out_snploc_stream $this_snp_id . ', ' . $snp_collection{$this_snp_id}{'pos'} . ', ' . $chr . "\n";
}
close $out_snploc_stream;



################
#subs
################
sub get_ancestral_allele{
	my ($info_string) = @_;
	my @info_fields = split(";", $info_string);
	if($info_fields[0] =~ /^AA=([A-Za-z]+)/){
		return uc($1);
	}else{
		print "Warning: ancestral allele info: $info_string not present/recognised.\n";
		return '-';
	}
}
#-------------------------------------------------------------------------------


#the desired output is a hash with key the sample name, and value a string as follows:
#XX
#?? if no data available
#example input
#T C 1/1:0.8888      0/1:0.8888      1/1:0.8888      0/1:0.8888      0/1:0.8004 
sub get_sample_to_genotype_hash{
	my ($my_vcf_ref, $my_vcf_alt, %my_sample_to_gt) = @_;
	#my ($my_vcf_ref, $my_vcf_alt, $my_vcf_anc, %my_sample_to_gt) = @_;
	my %sample_to_genotype;
	
#	my %allele_map;
#	#------------------------------
#	#first compare ref and ancestral
#	#------------------------------
#	#TODO very simplistic, needs more work
#	if( ($my_vcf_anc eq $my_vcf_ref) or ($my_vcf_anc eq '-') ){
#		%allele_map = (
#			'0' => $my_vcf_ref,
#			'1' => $my_vcf_alt 
#		);		
#	}elsif($my_vcf_anc eq $my_vcf_alt){
#		#then we need to look at the alternate
#		%allele_map = (
#			'0' => $my_vcf_alt,
#			'1' => $my_vcf_ref 
#		);	
#	}else{
#		print "ERROR: ancestral, reference and alternate differ: $my_vcf_anc, $my_vcf_ref, $my_vcf_alt. Aborting..\n";
#		exit -1;
#	}

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
			my $output = $allele_map{$1} . $allele_map{$2};
			$sample_to_genotype{$item} = $output;
		}else{
			#genotype no present
			$sample_to_genotype{$item} = 'NN';
		}
	}

	#my @temp;
	#foreach my $item (sort keys %sample_to_genotype){
	#	push(@temp, $sample_to_genotype{$item});
	#}
	#return join(", ", @temp);
	return \%sample_to_genotype;
}

#copied from do_diffbind_per_peak_analysis.pl
#maybe useful?
#sub get_sample_to_gtlabel_hash{
#	my ($my_vcf_ref, $my_vcf_alt_string, $my_vcf_anc, %my_sample_to_gt) = @_;
#	my %my_sample_to_gtlabels;
#	my %gt_to_label;
#	
#	#------------------------------
#	#first compare ref and ancestral
#	#------------------------------
#	if( ($my_vcf_anc eq $my_vcf_ref) or ($my_vcf_anc eq '-') ){
#		#then 0|0 gets a HOMREF label while all the rest gets a NON-HOMREF label
#		$gt_to_label{'0|0'} = 1; #all the rest won't be found in this hash and gets a NON-HOMREF
#		$gt_to_label{'0/0'} = 1;
#	}elsif($my_vcf_anc ne $my_vcf_ref){
#		#then we need to look at the alternate
#		#we first need to establish how many alternate genotypes are there
#		if($my_vcf_alt_string =~ /\w{1}/){ #only one------------------------------
#			if($my_vcf_alt_string eq $my_vcf_anc){
#				#then 1|1 gets a HOMREF label while all the rest gets a NON-HOMREF label
#				$gt_to_label{'1|1'} = 1; #all the rest won't be found in this hash and gets a NON-HOMREF
#				$gt_to_label{'1/1'} = 1;
#			}else{
#				print "Error: ref,alternate and ancestral differ: $my_vcf_ref, $my_vcf_alt_string, $my_vcf_anc. Skipping..\n";
#				return undef;
#			}
#		}elsif($my_vcf_alt_string =~ /(\w\W)+/){ #multiple alternate A,G--------
#			my $pos = 1; my $found;
#			my @alternate_gts = split(',', $my_vcf_alt_string);
#			foreach my $altg (@alternate_gts){
#				if($altg eq $my_vcf_anc){
#					#then $pos|$pos gets a HOMREF label, while all the rest gets a NON-HOMREF label
#					my $pos_pos = $pos . '|' . $pos;
#					$gt_to_label{$pos_pos} = 1; #all the rest won't be found in this hash and gets a NON-HOMREF
#					$pos_pos = $pos . '/' . $pos;
#					$gt_to_label{$pos_pos} = 1;
#					$found = 1;
#					last;
#				}else{
#					$pos += 1;
#				}
#			}
#			if(!$found){
#				print "Error: ancestral: $my_vcf_anc non found in list of alternates: $my_vcf_alt_string. Skipping..\n";
#				return undef;				
#			}
#		}else{
#			#should never get in here
#			print "Error: ref different from ancestral, and alternate string not recognised: $my_vcf_alt_string. Skipping..\n";
#			return undef;
#		}
#	}
#	#now we have a hash mapping genotype code to genotype label. We need to assign one of these labels for each sample.
#	foreach my $item (sort keys %my_sample_to_gt){
#		my @genotype_fields = split(':', $my_sample_to_gt{$item});
#		my $sample_gt = $genotype_fields[0];
#		#check that the genotype is actually present: 
#		unless($sample_gt =~ /\d{1}[\/|\|]\d{1}/){
#			print "Warning: genotype for sample: $item not recognised/not present: $sample_gt. Skipping sample..\n";
#			next;
#		}
#		#if all is correct, create label
#		if($gt_to_label{$sample_gt}){# then it has the homozygous ancestral
#			$my_sample_to_gtlabels{$item} = $LABEL_HOM;
#		}else{
#			$my_sample_to_gtlabels{$item} = $LABEL_NONHOM;
#		}
#	}
#	return \%my_sample_to_gtlabels;
#}