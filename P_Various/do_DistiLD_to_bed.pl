#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
#use String::Approx 'amatch';

#INFO
#THIS uses all annotation in distilLD as is, provided the term appears for at least #minsnp SNPs and the snp is genome-wide significant

#INFO 15/11
#This converts the DistiLD snps.tsv.gz table found here
#http://distild.jensenlab.org/download.html
#into an annotated BED file for usage with GAT

#INFO 15/10
#Added input argument to select minimum p value threshold. Default is 5x10-8
#only intervals passing this threshold will be included

#INPUT format
#0 PubMed ID of GWAS study
#1 Reference SNP (rs) number
#2 P-value
#3 Linkage disequilibrium (LD) block
#4 Ensembl genes in LD block
#5 Short description of GWAS study
#6 ICD10 codes (if applicable)

#eg
#0 21529783
#1 rs9556711
#2 8E-7
#3 chr13:97521524-98312811
#4 ENSG00000125249;ENSG00000125249;ENSG00000125249;ENSG00000165621;ENSG00000139793
#5 Alcoholism 12-month weekly alcohol consumption
#6 F10.2

#turn this into something of the form:
#ATTENTION! You need to do a mergebed after
#CHROM	START	END	NAME	SCORE
#3-1	3-2	3-3	1(5)	2 

#OUTPUT 1 files:
#segment file for GAT

#You need to create a temp bed file per phenotype
#use a merge bed on it to merge overlapping features
#then output the union bed file

#classify based on field 5
#notice that many times you will have things like "type I diabetes" and "type 1 diabetes". Equalise

#TODO
#my $PVAL_THRS = '5E-8'; #change if needed

my $infile;
my $min_snp_number;
my $phenotype;
#GetOptions(
#        'i=s'         =>\$infile,
#        'minsnp=i'    =>\$min_snp_number,
#        'ph=s'        =>\$phenotype
#);
GetOptions(
        'i=s'         =>\$infile,
        'ph=s'        =>\$phenotype
);

if(!$infile){
     print "USAGE: do_DistiLD_to_bed.pl -i=<INFILE> -minsnp=<MIN_SNP_NUMBER> -ph=<PHENOTYPE>\n";
     print "<INFILE> DistiLD TSV file\n";
     print "<MIN_SNP_NUMBER> minimum number of SNP per phenotype to consider for output\n";
     print "<PHENOTYPE> (optional) output contains data for the indicated phenotype only (default: all phenotypes, compatibly with <MIN_SNP_NUMBER>)\n";
     exit 1;
}
#if(!$min_snp_number){
#     print "USAGE: do_DistiLD_to_bed.pl -i=<INFILE> -minsnp=<MIN_SNP_NUMBER> -ph=<PHENOTYPE>\n";
#     print "<INFILE> DistiLD TSV file\n";
#     print "<MIN_SNP_NUMBER> minimum number of SNP per phenotype to consider for output\n";
#     print "<PHENOTYPE> (optional) output contains data for the indicated phenotype only (default: all phenotypes, compatibly with <MIN_SNP_NUMBER>)\n";
#     exit 1;
#}

my $TOOL_PATH       = '/net/isi-scratch/giuseppe/tools';
my $BEDTOOLS  = $TOOL_PATH . '/bedtools2-2.20.1/bin';

my ($filename, $directory) = fileparse($infile);
$filename =~ s/(.*)\..*/$1/;
my $outfile = $directory . "\/" .  $filename . '_segments.bed';

#1st pass
#build phenotype catalog and phenotype->snp-number map
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
my %phenotype_map;
while(<$instream>){
	chomp;
	my $description = (split /\t/)[5];
	#next if ($description =~ /quantitative trait/);
	next if ($description =~ /^response to/);
	#need to correct some ambiguities
	#type I diabetes -> type 1 diabetes
	#type II diabetes -> type 2 diabetes
	$description =~ s/type i/type 1/i   if(lc($description) =~ /type i diabetes/);
	$description =~ s/type ii/type 2/i  if(lc($description) =~ /type ii diabetes/);
	#alzheimers disease -> alzheimer's disease
	$description =~ s/alzheimers/alzheimer's/i   if(lc($description) =~ /alzheimers/);
	#parkinsons -> parkinson's
	$description =~ s/parkinsons/parkinson's/i   if(lc($description) =~ /parkinsons/);
	
	$description =~ s/body mass index bmi/body mass index/i   if(lc($description) =~ /body mass index bmi/);
	$description =~ s/sporadic amyotrophic lateral sclerosis als/sporadic amyotrophic lateral sclerosis/i   if(lc($description) =~ /sporadic amyotrophic lateral sclerosis als/);
	$description =~ s/systemic lupus erythematosus sle/systemic lupus erythematosus/i   if(lc($description) =~ /systemic lupus erythematosus sle/);
	$description =~ s/testicular germ cell tumor tgct/testicular germ cell tumor/i   if(lc($description) =~ /testicular germ cell tumor tgct/);
	$description =~ s/type 1 diabetes t1d/type 1 diabetes/i   if(lc($description) =~ /type 1 diabetes t1d/);
	$description =~ s/type 2 diabetes mellitus/type 2 diabetes/i   if(lc($description) =~ /type 2 diabetes mellitus/);
	$description =~ s/age-related macular degeneration amd/age-related macular degeneration/i   if(lc($description) =~ /age-related macular degeneration amd/);
	$description =~ s/amyotrophic lateral sclerosis als/amyotrophic lateral sclerosis/i   if(lc($description) =~ /amyotrophic lateral sclerosis als/);
	$description =~ s/autism spectrum disorders asds/autism spectrum disorders/i   if(lc($description) =~ /autism spectrum disorders asds/);

	if($phenotype){ #then only do it for the required phenotype
		if(lc($description) eq lc($phenotype)){ #here we probably need something more powerful
			if(!$phenotype_map{lc($description)}){
				$phenotype_map{lc($description)} = 1;
			}else{
				$phenotype_map{lc($description)} += 1;		
			}	
		}else{
			next;
		}
	}

	if(!$phenotype_map{lc($description)}){
		$phenotype_map{lc($description)} = 1;
	}else{
		$phenotype_map{lc($description)} += 1;		
	}
}
close $instream;

if($phenotype){
	if(!$phenotype_map{lc($phenotype)}){
		print "Unable to find phenotype $phenotype in DistiLD list. Aborting..\n";
		exit;
	}	
}

#create one bed for each of the terms above, if the term has > $min_snp_number
my %collective_bed;
foreach my $item (sort keys %phenotype_map){

	#next if($phenotype_map{$item} < $min_snp_number);
	
	my $outfile_single_phenotype           = $directory . "\/" .  $filename . '_phenotype_' . '_temp.bed';
	my $outfile_single_phenotype_processed = $directory . "\/" .  $filename . '_phenotype_' . '_temp_out.bed';
	open ($instream,  q{<}, $infile) or die("Unable to open $infile : $!");
	open (my $temp_stream,  q{>}, $outfile_single_phenotype) or die("Unable to open $outfile_single_phenotype : $!");
	while(<$instream>){
		chomp;
		next if($_ eq '');
		my $chr; my $start; my $end;
		my $ld_block;  #3
		my $desc;      #5

		($ld_block, $desc) = (split /\t/)[3,5];
		if(lc($desc) eq $item){
			#print $item . ": " . $description . "\n";
			if ($ld_block =~ /(chr\w+):(\d+)-(\d+)/){
				$chr = $1; $start = $2; $end = $3;
			}else{
				print STDERR "LD block: $ld_block in line $_ is in unrecognisable format.\n";
				#exit -1;
				next;
			}
			my $line = $chr . "\t" . $start . "\t" . $end  . "\t" . $item .  "\n";
			print $temp_stream 	$line;
		}
	}
	close $instream;
	close $temp_stream;
	#bedtools merge [OPTIONS] -i <BED/GFF/VCF>
	system "sort -k1,1V -k2,2g $outfile_single_phenotype | $BEDTOOLS/bedtools merge -i stdin > $outfile_single_phenotype_processed";
	unlink $outfile_single_phenotype;
	open ($temp_stream,  q{<}, $outfile_single_phenotype_processed) or die("Unable to open $outfile_single_phenotype_processed : $!");
	while(<$temp_stream>){
		chomp;
		next if($_ eq '');
		my $data_line = $_ . "\t" . $item . "\n";
		$collective_bed{$data_line} = 1;
	}
	close $temp_stream;
	unlink $outfile_single_phenotype_processed;
	#run bedtools, open file, save data in %collective_bed, delete temp file
}
open (my $segment_stream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
foreach my $item (sort keys %collective_bed){ print $segment_stream $item; }
close $segment_stream;


#open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
#
##save into hash to remove duplicates
#my %segments;
#my %descriptions;
#while(<$instream>){
#	chomp;
#	next if($_ eq '');
#	#bed fields
#	my $chr; my $start; my $end;
#	my $ld_block;  #3
#	my $snp_id;    #1
#	my $desc;      #5
#	my $p_val;     #2  
#
#	($snp_id, $p_val, $ld_block, $desc) = (split /\t/)[1,2,3,5];
#	if ($ld_block =~ /(chr\w+):(\d+)-(\d+)/){
#		$chr = $1; $start = $2; $end = $3;
#	}else{
#		print STDERR "LD block: $ld_block in line $_ is in unrecognisable format.\n";
#		#exit -1;
#		next;
#	}
#	#build name field
#	#my $name = $snp_id . ' (' . $desc . ')';
#	#my $line = $chr . "\t" . $start . "\t" . $end . "\t" . $name . "\t" . $p_val . "\n";
#	my $line = $chr . "\t" . $start . "\t" . $end . "\t" . $snp_id . "\t" . $p_val . "\n";
#	$segments{$line} = 1;
#	$descriptions{$snp_id} =  $desc;
#}
#close $instream;
#
#open (my $segment_stream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
#foreach my $item (keys %segments){ print $segment_stream $item; }
#close $segment_stream;
#
#open (my $descriptions_stream,  q{>}, $outfile_descriptions) or die("Unable to open $outfile_descriptions : $!");
#print $descriptions_stream "snp_id\tdescription\n";
#foreach my $item (keys %descriptions){ print $descriptions_stream $item . "\t" . $descriptions{$item} . "\n"; }
#close $descriptions_stream;



#string::approx
#sub _fuzzy_match {
#     my ($sig, $id_to_test) = @_;
# 
#     return amatch($sig, [ # this array sets match options:
#                              "i",    # match case-insensitively
#                         ], $id_to_test);
#}
