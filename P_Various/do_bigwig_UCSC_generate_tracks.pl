#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
#use PerlIO::gzip;

#THIS SCRIPT ASSUMES YOU HAVE GENERATED YOUR BIGWIG FILES USING MACS2 (bedgraph, do_bedgraph_to_bigwig.sh) and copied them to the public server at
#/home/giuseppe/public_html

#script to create UCSC bigwig track. 
#it is supposed to visualise data around a snp of interest
#samples will be sorted by genotype

#Cut and paste the output of this file into UCSC (custom tracks) to load them. You will probably need to reload them with an adjusted maximum height once you have #visualised a peak of interest.

#examples
#track name=NA06986  description="NA06986 - CEU, A/A" autoScale=off windowingFunction=mean+whiskers viewLimits=0:8 maxHeightPixels=100:100:5 smoothingWindow=5 bigDataUrl=http://wwwfgu.anat.ox.ac.uk/~giuseppe/NA06986.bw type=bigWig
#track type=bigWig name=proteinA smoothingWindow=4 color=123,100,50 autoScale=on viewLimits=1:200 visibility=full windowingFunction=maximum bigDataUrl=https://projects/files/file.bw

#Inputs
#sample names
#genotypes for a particular peak? (so that you can provide an ordering)
#position for a particular peak? (so that you can provide a section)
#http link

my $VIEWLIMITS;
my $SNP_rsID;
my $gt_choice;
my $GENOTYPE_FILE;

GetOptions(
        'snp=s'		   =>\$SNP_rsID,
		'gt=s'         =>\$gt_choice,
        'vl=f'         =>\$VIEWLIMITS
);
if(!$SNP_rsID){
     print "USAGE: do_ASSOCIATION_snptest.pl -snp=<SNP_RSID> -gt=<1kg|impute> -vl=<VIEW_LIMITS>\n";
     print "<SNP_RSID> id of the snp implicated in the association\n";
     print "<1kg|impute> whether to search on the 1kg vcf or the imputed 1kg vcf\n";
     print "<VIEW_LIMITS> number indicated the maximum view limit (eg 8)\n";
     exit 1;
}
if(!$VIEWLIMITS){
     print "USAGE: do_ASSOCIATION_snptest.pl -snp=<SNP_RSID> -gt=<1kg|impute> -vl=<VIEW_LIMITS>\n";
     print "<SNP_RSID> id of the snp implicated in the association\n";
     print "<1kg|impute> whether to search on the 1kg vcf or the imputed 1kg vcf\n";
     print "<VIEW_LIMITS> number indicated the maximum view limit (eg 8)\n";
     exit 1;
}

$VIEWLIMITS = '0:' . $VIEWLIMITS;
#change as needed
my $DATA_DIR = '/home/giuseppe/public_html';
my $SMOOTHING_WIN = 5;
my $AUTOSCALE = 'off';
my $VISIBILITY = 'full';
my $WINDOWING_F = 'mean+whiskers';
my $MAX_HEIGHT_P = '100:100:5';

#I want to pick the line with the snp here and get the genotypes
#I will need the snp rsID

if($gt_choice eq '1kg'){
	$GENOTYPE_FILE = "/net/isi-scratch/giuseppe/VDR/VARIANTS/1000g_hg19_b37/ALL.phase1_release_v3.20101123.snps_indels_svs.genotypes_FIFTEEN_1K_GEN_hg19.vcf";
}elsif($gt_choice eq 'impute'){
	$GENOTYPE_FILE = "/net/isi-scratch/giuseppe/VDR/VARIANTS/IMPUTATION/d_IMPUTE_out_chksize_4mb/IMPUTE_1kg_to_hapmap_autosomes_hg19.vcf";
}else{
	print "ERROR: -gt field: $gt_choice not recognised. Aborting..\n";
	exit -1;
}

print "READING bigwig data from directory: $DATA_DIR\n";
print "USING genotype file: $GENOTYPE_FILE\n\n"; 
##############
#Genotype processing
##############
#get header to map samples to gts (copied from do_association_BIMBAM.pl)
my @c_names;
my %sample_to_gt;
open (my $instream,     q{<}, $GENOTYPE_FILE) or die("Unable to open $GENOTYPE_FILE : $!");
while(<$instream>){
	chomp;
	my @cols;
	next if ($_ =~ /^\#\#/); #metadata	
	if($_ =~ /^\#CHROM/){ #header
			@cols = split("\t", $_);
		    @c_names = @cols[9..(@cols-1)];
		    last;
	}
}
close $instream;

#use zgrep
#my %sample_to_gt;
#my $search_string = '"^\#CHROM"';
#my $vcf_header_line = `zgrep -P $search_string $GENOTYPE_FILE`;
#my @cols = split("\t", $vcf_header_line);
#my @c_names = @cols[9..(@cols-1)];

#search snp of interest against genotypes
my $search_string = '"\W' . $SNP_rsID . '\W"';
my $vcf_line = `grep -P $search_string $GENOTYPE_FILE`;
#map sample names to genotypes
my @vcf_line = split(/\t/, $vcf_line);
my @genotypes = @vcf_line[9..$#vcf_line]; 
#map all samples in the vcf to their genotype:
#(for some of these I won't have normalised reads)
@sample_to_gt{@c_names} = @genotypes;

my $snp_chr = $vcf_line[0]; 
my $snp_loc = $vcf_line[1];
my $snp_ref = $vcf_line[3];
my $snp_alt = $vcf_line[4];

my $gt_hr  = $snp_ref . $snp_ref;
my $gt_het = $snp_ref . $snp_alt;
my $gt_hnr = $snp_alt . $snp_alt;

#the next function assumes there is ONLY ONE ALTERNATE ALLELE (not a list)
#if this is not the case, next
if($snp_alt =~ /(\w\W)+/){
	print "WARNING: alternate allele string: $snp_alt seem to contain multiple alternate alleles. Aborting..\n";
	exit;
}
my $s_to_g = get_sample_to_genotype_hash($snp_ref, $snp_alt, %sample_to_gt);


##order by genotype
###############
#get track name and description from filename
###############
chdir $DATA_DIR;
my @files = <*.bw>;

my %sample_homref_to_line;
my %sample_het_to_line;
my %sample_homnr_to_line;

foreach my $infile (@files) {
	my $NAME;
	my $DESCRIPTION;
	my $SERVER_URL = 'http://wwwfgu.anat.ox.ac.uk/~giuseppe/';

	if($infile =~ /(NA\d{5}).+/){
		$NAME = $1;
	}else{
		print "Attention: big wig file name: $infile not recognised. Aborting\n";
		exit -1;
	}
	$SERVER_URL = $SERVER_URL . $infile;

	# common part of the line - description changes
	my $line =  'track type=bigWig' .
				' name='        . $NAME          .   
				' autoScale='         . $AUTOSCALE     .   
				' windowingFunction=' . $WINDOWING_F   .
				' visibility=full'                     .
				' gridDefault=off'                     . 
				' viewLimits='        . $VIEWLIMITS    . 
				' maxHeightPixels='   . $MAX_HEIGHT_P  . 
				' smoothingWindow='   . $SMOOTHING_WIN .
				' bigDataUrl='        . $SERVER_URL;

	#check the genotype for this sample
	my $gt = $$s_to_g{$NAME};
	if($gt eq $gt_hr){
		$DESCRIPTION = $NAME . ' - ' . $SNP_rsID . ' - ' . $snp_chr . '-' . $snp_loc . ' '  .  $gt . ' - ' .  'HOM-REF';
		$line =  $line . ' description="'      . $DESCRIPTION   . '"'; 
		$sample_homref_to_line{$NAME} = $line;
	}elsif($gt eq $gt_het){
		$DESCRIPTION = $NAME . ' - ' . $SNP_rsID . ' - ' . $snp_chr . '-' . $snp_loc . ' '  .  $gt . ' - ' .  'HET';
		$line =  $line . ' description="'      . $DESCRIPTION   . '"'; 
		$sample_het_to_line{$NAME} = $line;
	}elsif($gt eq $gt_hnr){
		$DESCRIPTION = $NAME . ' - ' . $SNP_rsID . ' - ' . $snp_chr . '-' . $snp_loc . ' '  .  $gt . ' - ' .  'HOM-NR';
		$line =  $line . ' description="'      . $DESCRIPTION   . '"'; 
		$sample_homnr_to_line{$NAME} = $line;
	}elsif($gt eq 'NN'){
		print "Warning: sample $NAME has genotype NN for this SNP. Skipping sample..\n";
		next;
	}
	else{
		print "Error: sample genotype: $gt does not match any allele configuration for this snp. Aborting.\n";
		exit;
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

#add a browser position track
#browser position chr21:33,031,597-33,041,570
my $browser_pos = $snp_chr . ':' . ($snp_loc - 300) . '-' . ($snp_loc + 300);
print 'browser position ' . $browser_pos . "\n";


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

