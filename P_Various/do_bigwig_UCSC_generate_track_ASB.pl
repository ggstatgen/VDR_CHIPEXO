#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
#use PerlIO::gzip;

#THIS SCRIPT ASSUMES YOU HAVE GENERATED YOUR BIGWIG FILES USING MACS2 (bedgraph, do_bedgraph_to_bigwig.sh) and copied them to the public server at
#/home/giuseppe/public_html

#script to create a UCSC bigwig track for looking at SNPs in significantly ASB events 


#Cut and paste the output of this file into UCSC (custom track) to load it. You will probably need to reload them with an adjusted maximum height once you have #visualised a peak of interest.


#track type=bigWig name=NA06986 autoScale=off windowingFunction=mean+whiskers visibility=full gridDefault=off viewLimits=0:5 maxHeightPixels=100:100:5 smoothingWindow=2 bigDataUrl=http://wwwfgu.anat.ox.ac.uk/~giuseppe/NA06986_d40_dupall.bw description="NA06986 - rs7094595 - 10-64415413"
#browser position 10:64415113-64415713



#examples
#track name=NA06986  description="NA06986 - CEU, A/A" autoScale=off windowingFunction=mean+whiskers viewLimits=0:8 maxHeightPixels=100:100:5 smoothingWindow=5 bigDataUrl=http://wwwfgu.anat.ox.ac.uk/~giuseppe/NA06986.bw type=bigWig
#track type=bigWig name=proteinA smoothingWindow=4 color=123,100,50 autoScale=on viewLimits=1:200 visibility=full windowingFunction=maximum bigDataUrl=https://projects/files/file.bw

#$GENOTYPE_FILE = "/net/isi-scratch/giuseppe/VDR/VARIANTS/1000g_hg19_b37/ALL.phase1_release_v3.20101123.snps_indels_svs.genotypes_FIFTEEN_1K_GEN_hg19.vcf";
#$GENOTYPE_FILE = "/net/isi-scratch/giuseppe/VDR/VARIANTS/IMPUTATION/d_IMPUTE_out_chksize_4mb/IMPUTE_1kg_to_hapmap_autosomes_hg19.vcf";

#Inputs
#sample name
#genotypes for a particular peak? (so that you can provide an ordering)
#position for a particular peak? (so that you can provide a section)
#http link

my $INPUT_VARIANTS_BASE;
my $PREFIX;
my $POSTFIX = '.vcf.gz';

my $VIEWLIMITS;
my $SNP_rsID;
my $sample_ID;
my $gt_choice;
my $GENOTYPE_FILE;

GetOptions(
        'snp=s'		   =>\$SNP_rsID,
        'sampleid=s'   =>\$sample_ID,
		'gt=s'         =>\$gt_choice,
        'vl=f'         =>\$VIEWLIMITS
);
if(!$SNP_rsID){
     print "USAGE: do_ASSOCIATION_snptest.pl -snp=<SNP_RSID> -sampleid=<ID> -gt=<1kg|impute> -vl=<VIEW_LIMITS>\n";
     print "<SNP_RSID> id of the snp implicated in the association\n";
     print "<ID> sample ID in NAxxxxx format \n";
     print "<1kg|impute> whether to search on the 1kg vcf or the imputed 1kg vcf\n";
     print "<VIEW_LIMITS> number indicated the maximum view limit (eg 8)\n";
     exit 1;
}
if(!$sample_ID){
     print "USAGE: do_ASSOCIATION_snptest.pl -snp=<SNP_RSID> -sampleid=<ID> -gt=<1kg|impute> -vl=<VIEW_LIMITS>\n";
     print "<SNP_RSID> id of the snp implicated in the association\n";
     print "<ID> sample ID in NAxxxxx format \n";
     print "<1kg|impute> whether to search on the 1kg vcf or the imputed 1kg vcf\n";
     print "<VIEW_LIMITS> number indicated the maximum view limit (eg 8)\n";
     exit 1;
}
if(!$gt_choice){
     print "USAGE: do_ASSOCIATION_snptest.pl -snp=<SNP_RSID> -sampleid=<ID> -gt=<1kg|impute> -vl=<VIEW_LIMITS>\n";
     print "<SNP_RSID> id of the snp implicated in the association\n";
     print "<ID> sample ID in NAxxxxx format \n";
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
my $SMOOTHING_WIN = 2;
my $AUTOSCALE = 'off';
my $VISIBILITY = 'full';
my $WINDOWING_F = 'mean+whiskers';
my $MAX_HEIGHT_P = '100:100:5';

#I want to pick the line with the snp here and get the genotypes
#I will need the snp rsID

if($gt_choice eq '1kg'){
	$INPUT_VARIANTS_BASE = '/net/isi-scratch/giuseppe/VDR/VARIANTS/1000g_hg19_b37/d_SPLIT_SINGLE_b37/';
	$PREFIX = 'ALL.phase1_release_v3.20101123.snps_indels_svs.genotypes_FIFTEEN_1K_GEN_';
}elsif($gt_choice eq 'impute'){
	$INPUT_VARIANTS_BASE = '/net/isi-scratch/giuseppe/VDR/VARIANTS/IMPUTATION/d_IMPUTE_out_chksize_4mb/d_ALLELESEQ/';
	$PREFIX = 'IMPUTE_1kg_to_hapmap_autosomes_g1k_v37_';
}else{
	print "ERROR: -gt field: $gt_choice not recognised. Aborting..\n";
	exit -1;
}
$GENOTYPE_FILE = $INPUT_VARIANTS_BASE . $PREFIX . $sample_ID . $POSTFIX;
	
print "READING bigwig data from directory: $DATA_DIR\n";
print "USING genotype file: $GENOTYPE_FILE\n\n"; 


#search snp of interest against genotypes
my $search_string = '"\W' . $SNP_rsID . '\W"';
my $vcf_line = `zgrep -P $search_string $GENOTYPE_FILE`;
#map sample names to genotypes
my @vcf_line = split(/\t/, $vcf_line); 

my $snp_chr = $vcf_line[0]; 
my $snp_loc = $vcf_line[1];
my $snp_ref = $vcf_line[3];
my $snp_alt = $vcf_line[4];


##order by genotype
###############
#get track name and description from filename
###############
chdir $DATA_DIR;
my @files = <*.bw>;

my %sample_to_line;
foreach my $infile (@files) {
	my $NAME;
	my $DESCRIPTION;
	my $SERVER_URL = 'http://wwwfgu.anat.ox.ac.uk/~giuseppe/';

	if($infile =~ /(NA\d{5}).+/){
		next if($1 ne $sample_ID);
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
		$DESCRIPTION = $NAME . ' - ' . $SNP_rsID . ' - ' . $snp_chr . '-' . $snp_loc;
		$line =  $line . ' description="'      . $DESCRIPTION   . '"'; 
		$sample_to_line{$NAME} = $line;
}
print $sample_to_line{$sample_ID}, "\n";

#add a browser position track
#browser position chr21:33,031,597-33,041,570
my $browser_pos = 'chr' . $snp_chr . ':' . ($snp_loc - 300) . '-' . ($snp_loc + 300);
print 'browser position ' . $browser_pos . "\n";

