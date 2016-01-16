#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
#use Statistics::Regression;

#old approach before testing regression based ones

#INFO 13/12
#This script carries out a per-peak Diffbind analysis of the samples
#Having ascertained that an INDIVIDUAL vs INDIVIDUAL analysis is not possible due to the absence of biological replicates,
#here I try to carry out an analysis at the peak level: the contrast is no more individual VS individual, but, given a peak, reference genotype VS non-reference genotype
#EDIT 13/12 - Chris suggests ancestral genotype vs non-ancestral genotype
#For each peak/SNP combination, can I split the samples in
#a) those with the reference genotype for the SNP (0/0) (REF)
#b) those with an alternative genotype for the SNP (0/1 or 1/1) (HET or HOM-NON-REF)

#INPUTS:
#1 peak file
#NOTE if you use the output from broadpeaks MACS2, you need to obtain the bed file from the .gappedPeak file using -cut f 1,2,3,4,5. The normal bed produced by macs2 lacks the 5th col.
#2 bam files for the samples
#3 disease snp file
#4 vcf with hapmap snps and genotypes

#intersect 1 and 2: this gives you a peak list with intersecting snps
#look up these snps in 3

#OUTPUT
#maybe a basic table telling which, if any peaks are differentially enriched? Or maybe a group of .csv tables to pass on to diffbind, one per peak

#ancestral alleles are defined as follows
#ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/technical/reference/ancestral_alignments/README

my $in_peak_bed;
my $in_disease_snp_bed;
my $in_bam;
my $in_genotypes_vcf;
my $datadir;
my $TOOL_PATH = '/net/isi-scratch/giuseppe/tools';
my $BEDTOOLS  = $TOOL_PATH . '/bedtools-2.17.0/bin';
#for now I hard code it
#my $GENOTYPES = '/net/isi-scratch/giuseppe/VDR/VARIANTS/1000g_hg19_b37/ALL.phase1_release_v3.20101123.snps_indels_svs.genotypes_FIFTEEN_1K_GEN_hg19.vcf'; #1000g
my $GENOTYPES = '/net/isi-scratch/giuseppe/VDR/VARIANTS/hapmap_liftover_hg19_b37/genotypes_CEUYRI_30_r27_nr.hg19_fwd.vcf'; #hapmap

my $LABEL_HOM    = 'HOM_ANCESTRAL';
my $LABEL_NONHOM = 'NON_HOM_ANCESTRAL';

#for the 1000g genotype vcf (4.1), format is GT:DS:GL
#GT : genotype, encoded as allele values separated by either of "/" or "|". 
#The allele values are 0 for the reference allele (what is in the REF field), 
#1 for the first allele listed in ALT, 
#2 for the second allele list in ALT and so on. 
#For diploid calls examples could be 0/1, 1|0, or 1/2, etc. 
#For haploid calls, e.g. on Y, male non-pseudoautosomal X, or mitochondrion, only one allele value should be given; a triploid call might look like 0/0/1. 
#If a call cannot be made for a sample at a given locus, "." should be specified for each missing allele in the GT field 
#(for example "./." for a diploid genotype and "." for haploid genotype). The meanings of the separators are as follows 
#(see the PS field below for more details on incorporating phasing information into the genotypes):
#
#   / : genotype unphased
#    | : genotype phased

#GL : genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible genotypes 
#given the set of alleles defined in the REF and ALT fields. In presence of the GT field the same ploidy is expected and
#the canonical order is used; without GT field, diploidy is assumed. If A is the allele in REF and B,C,... are the alleles
#as ordered in ALT, the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j.  In other words, for
#biallelic sites the ordering is: AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc.  For example: GT:GL 0/1:-323.03,-99.29,-802.53 (Floats)

#this was found in the GATK info http://gatkforums.broadinstitute.org/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk
#DS 	Were any of the samples downsampled because of too much coverage?

GetOptions(
        'p=s'     =>\$in_peak_bed,
        'bam=s'   =>\$in_bam,
        'dsnp=s'  =>\$in_disease_snp_bed,
        #'g=s'     =>\$in_genotypes_vcf  
);
if(!$in_peak_bed){
     print "USAGE: do_diffbind_per_peak_analysis.pl -p=<PEAK_FILE> -bam=<BAM_DIR> -dsnp=<DISEASE_SNP_FILE>\n";
     print "<PEAK_FILE> bed/narrowpeak file with total peaks to consider for differential analysis\n";
     print '<BAM_DIR> directory containing bam files for the samples for the analysis (in the form "VDR_SAMPLENAME_reads.bam")\n';
     print "<DISEASE_SNP_FILE> bed file containing the location of the GWAS phenotypes of interest\n";
     exit 1;
}
if(!$in_bam){
     print "USAGE: do_diffbind_per_peak_analysis.pl -p=<PEAK_FILE> -bam=<BAM_DIR> -dsnp=<DISEASE_SNP_FILE>\n";
     print "<PEAK_FILE> bed/narrowpeak file with total peaks to consider for differential analysis\n";
     print '<BAM_DIR> directory containing bam files for the samples for the analysis (in the form "VDR_SAMPLENAME_reads.bam")\n';
     print "<DISEASE_SNP_FILE> bed file containing the location of the GWAS phenotypes of interest\n";
     exit 1;
}
if(!$in_disease_snp_bed){
     print "USAGE: do_diffbind_per_peak_analysis.pl -p=<PEAK_FILE> -bam=<BAM_DIR> -dsnp=<DISEASE_SNP_FILE>\n";
     print "<PEAK_FILE> bed/narrowpeak file with total peaks to consider for differential analysis\n";
     print '<BAM_DIR> directory containing bam files for the samples for the analysis (in the form "VDR_SAMPLENAME_reads.bam")\n';
     print "<DISEASE_SNP_FILE> bed file containing the location of the GWAS phenotypes of interest\n";
     exit 1;
}
#if(!$in_genotypes_vcf){
#     print "USAGE: do_diffbind_per_peak_analysis.pl -p=<PEAK_FILE> -dsnp=<DISEASE_SNP_FILE> -g=<HAPMAP_VCF_FILE>\n";
#     print "<PEAK_FILE> bed/narrowpeak file with total peaks to consider for differential analysis\n";
#     print "<DISEASE_SNP_FILE> bed file containing the location of the GWAS phenotypes of interest\n";
#     print "<HAPMAP_VCF_FILE> vcf file containing genotypes for all the samples to carry out the diff analysis for <INFILE1>\n";
#     exit 1;
#}
#I want all the outputs in a directory whose format will be $PATH/filename_diffbind_by_peak
my($basename, $directory) = fileparse($in_peak_bed);
$basename =~ s/(.*)\..*/$1/;
if($directory eq "\.\/"){ 
	$datadir = 'd_diffbind_by_peak' . $basename;
}else{
	$datadir =  $directory . 'd_diffbind_by_peak_' .  $basename;
}
system "mkdir $datadir";

my $intervals_with_SNPs_file = $datadir . "\/" . $basename .  '_peak_snps.bed';
#BEDTOOLS:
#use bedtools intersect to annotate each peak with a snp with snp + genotype info
#wo option: 
#Write the original A and B entries plus the number of base
#		pairs of overlap between the two features.
#		- Overlaps restricted by -f and -r.
#		  Only A features with overlap are reported.
system "$BEDTOOLS/intersectBed -wo -a $in_peak_bed -b $in_disease_snp_bed  > $intervals_with_SNPs_file";
#output is like follows for macs2 BROAD:
#0 chr7
#1 37382185
#2 37382561
#3 /net/isi-scratch/giuseppe/VDR/POOLED_4/02_BAM_BOWTIE/d_CEU_YRI/MACS2_picard_merged.sorted_peak_420458
#4 654
#5 chr7
#6 37382464
#7 37382465
#8 rs60600003(Multiple Sclerosis)
#9 1
#output is like follows for gps/gem:
#0 chr1
#1 117100855
#2 117101055
#3 1:117100955
#4 267.0
#5 chr1
#6 117100956
#7 117100957
#8 rs1335532(Multiple sclerosis)
#9 1

#IF YOU USE MACS2 BROAD, THERE IS NO FIELD 4 (score)

#intervals with SNPs is now the main file
#it has one line per peak/snp combo.
#for each iteration, you want to create a csv file and run a diffbind analysis
#each cvs file will contain as many row as samples, and in each the DBA_CONDITION will be either HOMREF or NON-HOMREF
#Therefore I am running one analysis per (peak,snp) pair
#Need to correct for multiple conditions if any result comes out
open (my $instream,  q{<}, $intervals_with_SNPs_file) or die("Unable to open $intervals_with_SNPs_file : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	my %sample_to_gt;
	my $n_rep_hom = 0; 
	my $n_rep_non_hom = 0;
	
	#I will need the peak coordinates, peak signal and the snps (anything else?)
	#7 needs to be cleaned
	my ($chr, $start, $end, $name, $score, $snp_ph) = (split /\t/)[0,1,2,3,4,8];
	my $snp;
	if($snp_ph =~ /^(rs\d+)\(.+\)/){
		$snp = $1;
	}else{
		print "Disease-snp string not recognised. Skipping\n";
		next;
	}
	
	print "SNP $snp - creating peak data...\t";
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#1 create a temp bed file with only this peak, and an input file for diffbind
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#we need a short identifier for this peak and snp: chr_start_stop_snp?
	my $peak_id = $chr . '_' . $start . '_' . $end . '_' . $snp;
	my $outfile_singlepeak_bed      = $datadir . "\/" .  'peak_' .  $peak_id . '.bed';
	my $outfile_singlepeak_diffbind = $datadir . "\/" .  'peak_' .  $peak_id . '.csv';
	open (my $peak_stream,  q{>}, $outfile_singlepeak_bed) or die("Unable to open $outfile_singlepeak_bed : $!");
	$name =~ s/.+(peak_\d+)$/$1/; #/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_CEU_YRI_1kg16/MACS2_picard_merged.sorted_peak_1858279 > peak_1858279
	print $peak_stream $chr . "\t" . $start . "\t" . $end . "\t" . $name . "\t" . $score . "\n";
	close $peak_stream;
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#2 get header and cell names from vcf file
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	my ($columns, $cell_names) = get_vcf_header();
	
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#3 search for the snp in the vcf data and split the cells by genotype
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	print "Searching for snp in database...\t";
	my $snp_line = `grep -w -i $snp $GENOTYPES`;
	if(!$snp_line){
		print "**SNP $snp not found in the 1000g dataset**\n";
		next;
	}
	chomp $snp_line;
	my @snp_line = split(/\t/, $snp_line);
	my $gt_snp_id = $snp_line[2];
	unless($gt_snp_id =~ /^rs\d+/i){ #to check that something funny is not happening with the output from grep
		print "Unrecognised snp ID: $gt_snp_id. Skipping...\n";
		next;
	}
	my $vcf_ref = $snp_line[3];
	my $vcf_alt_string = $snp_line[4];
	my $vcf_info = $snp_line[7];
	my $vcf_ancestral = get_ancestral_allele($vcf_info);
	
	my @genotypes = @snp_line[9..(@$columns-1)]; 
	@sample_to_gt{@$cell_names} = @genotypes;
	
	print "Getting sample genotypes...\t";
	my $sample_to_gtlabels = get_sample_to_gtlabel_hash($vcf_ref, $vcf_alt_string, $vcf_ancestral, %sample_to_gt);
	if(!$sample_to_gtlabels){
		print "get_sample_to_gtlabel_hash(): unable to return a map for snp $snp. Skipping..\n";
		next;
	}
	#count number of items in each class (number of replicates)
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#4 build .csv file with an entry for each sample
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------
	#you will need:
	#1 sample names (from @cell_names)
	#2 Tissue = CEU/YRI - extract from sample names
	#3 Factor = 'VDR'
	#4 Condition = $LABEL_HOM / $LABEL_NONHOM
	#5 Replicate = 1..$n_rep_hom / 1..$n_rep_non_hom
	#6 bam file for this cell
	#7 bed file for this peak
	#sample data:
	#SampleID,Tissue,Factor,Condition,Replicate,bamReads,bamControl,Peaks,PeakFormat,ScoreCol,LowerBetter
	#GM06986,CEU,VDR,HOM_ANCESTRAL,01,POOLED/04_BAM_STAMPY_MINMAPQ_11/VDR_GM06986_reads.bam,,POOLED/04_BAM_STAMPY_MINMAPQ_11/d_gem_4_30/d_beds/GM06986_2_GEM_events_pp.bed,raw,5,F
	open (my $diffbind_stream,  q{>}, $outfile_singlepeak_diffbind) or die("Unable to open $outfile_singlepeak_diffbind : $!");
	#header
	print $diffbind_stream "SampleID,Tissue,Factor,Condition,Replicate,bamReads,bamControl,Peaks,PeakFormat,ScoreCol,LowerBetter\n";
	foreach my $item (sort keys %{$sample_to_gtlabels}){
		my $item_ethnicity; 
		my $diffbind_entry;
		#create bam file path (warning, very dataset-dependent)
		my $na_sample = $item;
		$item =~ s/NA(\d+)/GM$1/;
		my $bam_file = $in_bam . '/' . 'VDR_' . $item . '_reads.bam';
		#get ethnicity (warning, very dataset-dependent)
		if($item =~ /^GM19\d+/){
			$item_ethnicity = 'YRI';
		}else{
			$item_ethnicity = 'CEU';
		}
		
		if($$sample_to_gtlabels{$na_sample} eq $LABEL_HOM){
			$n_rep_hom++;
			$diffbind_entry = $item . ',' . $item_ethnicity . ',' . 'VDR' . ',' .  $LABEL_HOM . ',' .  $n_rep_hom     . ',' . $bam_file . ',,' . $outfile_singlepeak_bed . ',raw,5,F'; 			
		}
		if($$sample_to_gtlabels{$na_sample} eq $LABEL_NONHOM){
			$n_rep_non_hom++;
			$diffbind_entry = $item . ',' . $item_ethnicity . ',' . 'VDR' . ',' .  $LABEL_NONHOM . ',' . $n_rep_non_hom . ',' . $bam_file . ',,' . $outfile_singlepeak_bed . ',raw,5,F';
		}
		print $diffbind_stream $diffbind_entry, "\n";
	}
	close $diffbind_stream;
	print "Done.\n";
}
close $instream;


#=====
#subs
#=====
#-------------------------------------------------------------------------------
sub get_vcf_header{
	my @cols; my @c_names;
	
	open (my $instream_vcf,  q{<}, $GENOTYPES) or die("Unable to open $GENOTYPES : $!");
	while(<$instream_vcf>){
		chomp;
		next if($_ eq '');        #empty
		next if ($_ =~ /^##/);    #comment
		if($_ =~ /^#CHROM/){      #header
			#count the number of samples, get the samples' names, and move on
			@cols = split("\t", $_);
		    @c_names = @cols[9..(@cols-1)];
		    my $cell_number = @cols - 9;
			#print "vcf contains genotypes for $cell_number samples.\n"; 
			last;
		}
	}
	close $instream_vcf;
	return (\@cols, \@c_names);
}
#-------------------------------------------------------------------------------
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
#when available, I will consider the ancestral allele as my reference allele. Therefore I will call "HOM_REF" genotypes where both alleles have the ancestral nucleotide
#purpose - for each cell, this needs to read the first field of the genotype info 0|1:1.000:-3.30,-0.00,-4.22 and assign a 
#label to the sample based on the ancestral (if present) or reference (if ancestral not present)
#example 20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
#vcf-ref = A, vcf-alt = G,T, vcf-anc=T, gt1 = 1|2, gt2 = 2|1
#Therefore:
#case 1: vcf-anc == vcf-ref -> all gt = 0|0 are HOMREF (HOMANC) and all gt != 0|0 are NON-HOMREF (NON-HOMANC)
#case 2: vcf-anc != vcf-ref. Then we only care about the vcf-anc. Check vcf-alt. 
#is vcf-alt only one char? YES: compare to vcf-anc. Are they the same? YES -> all gt = 1|1 are HOMANC and all gt != 0|0 are NON-HOMANC
sub get_sample_to_gtlabel_hash{
	my ($my_vcf_ref, $my_vcf_alt_string, $my_vcf_anc, %my_sample_to_gt) = @_;
	my %my_sample_to_gtlabels;
	my %gt_to_label;
	
	#------------------------------
	#first compare ref and ancestral
	#------------------------------
	if( ($my_vcf_anc eq $my_vcf_ref) or ($my_vcf_anc eq '-') ){
		#then 0|0 gets a HOMREF label while all the rest gets a NON-HOMREF label
		$gt_to_label{'0|0'} = 1; #all the rest won't be found in this hash and gets a NON-HOMREF
		$gt_to_label{'0/0'} = 1;
	}elsif($my_vcf_anc ne $my_vcf_ref){
		#then we need to look at the alternate
		#we first need to establish how many alternate genotypes are there
		if($my_vcf_alt_string =~ /\w{1}/){ #only one------------------------------
			if($my_vcf_alt_string eq $my_vcf_anc){
				#then 1|1 gets a HOMREF label while all the rest gets a NON-HOMREF label
				$gt_to_label{'1|1'} = 1; #all the rest won't be found in this hash and gets a NON-HOMREF
				$gt_to_label{'1/1'} = 1;
			}else{
				print "Error: ref,alternate and ancestral differ: $my_vcf_ref, $my_vcf_alt_string, $my_vcf_anc. Skipping..\n";
				return undef;
			}
		}elsif($my_vcf_alt_string =~ /(\w\W)+/){ #multiple alternate A,G--------
			my $pos = 1; my $found;
			my @alternate_gts = split(',', $my_vcf_alt_string);
			foreach my $altg (@alternate_gts){
				if($altg eq $my_vcf_anc){
					#then $pos|$pos gets a HOMREF label, while all the rest gets a NON-HOMREF label
					my $pos_pos = $pos . '|' . $pos;
					$gt_to_label{$pos_pos} = 1; #all the rest won't be found in this hash and gets a NON-HOMREF
					$pos_pos = $pos . '/' . $pos;
					$gt_to_label{$pos_pos} = 1;
					$found = 1;
					last;
				}else{
					$pos += 1;
				}
			}
			if(!$found){
				print "Error: ancestral: $my_vcf_anc non found in list of alternates: $my_vcf_alt_string. Skipping..\n";
				return undef;				
			}
		}else{
			#should never get in here
			print "Error: ref different from ancestral, and alternate string not recognised: $my_vcf_alt_string. Skipping..\n";
			return undef;
		}
	}
	#now we have a hash mapping genotype code to genotype label. We need to assign one of these labels for each sample.
	foreach my $item (sort keys %my_sample_to_gt){
		my @genotype_fields = split(':', $my_sample_to_gt{$item});
		my $sample_gt = $genotype_fields[0];
		#check that the genotype is actually present: 
		unless($sample_gt =~ /\d{1}[\/|\|]\d{1}/){
			print "Warning: genotype for sample: $item not recognised/not present: $sample_gt. Skipping sample..\n";
			next;
		}
		#if all is correct, create label
		if($gt_to_label{$sample_gt}){# then it has the homozygous ancestral
			$my_sample_to_gtlabels{$item} = $LABEL_HOM;
		}else{
			$my_sample_to_gtlabels{$item} = $LABEL_NONHOM;
		}
	}
	return \%my_sample_to_gtlabels;
}