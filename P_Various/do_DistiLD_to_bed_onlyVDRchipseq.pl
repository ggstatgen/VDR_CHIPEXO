#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
#use String::Approx 'amatch';

#INFO
#THIS "clusters" similar terms based on a dictionary search with "reference phenotypes" found on the VDR chip-seq paper.
#I simplify the distilLD annotation to group together all terms eg containing "multiple sclerosis"


#INFO 15/11
#This converts the DistiLD .txt table found here
#http://distild.jensenlab.org/download.html
#into an annotated BED file for usage with GAT

#Andreas confirmed today that they don't use anymore the GWAS catalog (extending the snp position by 150KB left and right and turning it into a bed) but use this instead
#because it contains more precise LD blocks

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

#OUTPUT 2 files:
#segment file for GAT --segments
#annotation name for GAT --descriptions

#You need to create a temp bed file per phenotype
#use a merge bed on it to merge overlapping features
#then output the union bed file


#classify based on field 5

my $infile;
GetOptions(
        'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_DistiLD_to_bed.pl -i=<INFILE>\n";
     print "<INFILE> DistiLD TSV file\n";
     exit 1;
}
my $TOOL_PATH       = '/net/isi-scratch/giuseppe/tools';
my $BEDTOOLS  = $TOOL_PATH . '/bedtools2-2.19.0/bin';

my ($filename, $directory) = fileparse($infile);
$filename =~ s/(.*)\..*/$1/;

my $outfile_phenotype_list = $directory . "\/" .  $filename . '_phenotype_list.txt';
my $outfile_segments = $directory . "\/" .  $filename . '_segments.bed';

#this is the list of phenotypes used by andreas in the VDR genome res paper
#http://genome.cshlp.org/content/20/10/1352.full#methods
my @reference_phenotypes = (
"acute myeloid leukemia",
"body mass index",
"aging",
"hiv",
"alcohol dependence",
"alzheimer's disease",
"amyotrophic lateral sclerosis",
"ankylosing spondylitis",
"asthma",
"atrial fibrillation",
"attention deficit hyperactivity disorder",
"bipolar disorder",
"blood pressure",
"bone mineral density",
"breast cancer",
"celiac disease",
"cholesterol",
"chronic lymphocytic leukemia",
"colorectal cancer",
"coronary artery disease",
"crohn's disease",
"fasting glucose",
"hair color",
"hdl cholesterol",
"height",
"ldl cholesterol",
"leprosy",
"lung cancer",
"melanoma",
"multiple sclerosis",
"myocardial infarction",
"pancreatic cancer",
"Parkinson's disease",
"primary biliary cirrhosis",
"\^prostate cancer",
"psoriasis",
"qt interval",
"restless legs syndrome",
"\^rheumatoid arthritis",
"schizophrenia",
"stroke",
"systemic lupus erythematosus",
"tanning",
"\^triglycerides",
"type 2 diabetes",
"ulcerative colitis"
);

my @reference_phenotypes_DistiLD = (
"acute lymphoblastic leukemia",
"aging\$", #
"aids",
"alcohol dependence",
"alzheimer",
"amyotrophic lateral sclerosis",
"ankylosing spondylitis",
"asthma\$",
"atrial fibrillation",
"attention deficit hyperactivity disorder",
"bipolar disorder",
"black vs. blond hair color",
"blood pressure",
"body mass index",
"bone mineral density",
"breast cancer",
"celiac disease",
"\^cholesterol",
"chronic lymphocytic leukemia",
"colorectal cancer",
"coronary heart disease",
"crohn's disease",
"fasting plasma glucose",
"hair",
"hdl cholesterol",
"height",
"hiv-1",
"iris",
"ldl cholesterol",
"leprosy",
"lung cancer",
"melanoma",
"multiple sclerosis",
"myocardial infarction",
"pancreatic cancer",
"parkinson",
"primary biliary cirrhosis",
"prostate cancer",
"prostate cancer",
"psoriasis",
"psoriatic arthritis",
"qt interval",
"rheumatoid arthritis",
"schizophrenia",
"stroke",
"systemic lupus erythematosus",
"tanning",
"\^triglycerides",
#"diabetes",
"type 1 diabetes",
"type 2 diabetes",
"type i diabetes",
"type ii diabetes",
"\^ulcerative colitis",
);

#create one bed for each of the terms above
my %collective_bed;
foreach my $item (@reference_phenotypes_DistiLD){
	my $outfile_single_phenotype           = $directory . "\/" .  $filename . '_phenotype_' . '_temp.bed';
	my $outfile_single_phenotype_processed = $directory . "\/" .  $filename . '_phenotype_' . '_temp_out.bed';
	open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
	open (my $temp_stream,  q{>}, $outfile_single_phenotype) or die("Unable to open $outfile_single_phenotype : $!");
	while(<$instream>){
		chomp;
		next if($_ eq '');
		my $chr; my $start; my $end;
		my $ld_block;  #3
		my $desc;      #5

		($ld_block, $desc) = (split /\t/)[3,5];
		$desc =~ s/type i/type 1/i   if(lc($desc) =~ /type i diabetes/);
		$desc =~ s/type ii/type 2/i  if(lc($desc) =~ /type ii diabetes/);
		#alzheimers disease -> alzheimer's disease
		$desc =~ s/alzheimers/alzheimer's/i   if(lc($desc) =~ /alzheimers/);
		#parkinsons -> parkinson's
		$desc =~ s/parkinsons/parkinson's/i   if(lc($desc) =~ /parkinsons/);
		
		if(lc($desc) =~ /$item/){
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
	open (my $temp_stream_processed,  q{<}, $outfile_single_phenotype_processed) or die("Unable to open $outfile_single_phenotype_processed : $!");
	while(<$temp_stream_processed>){
		chomp;
		next if($_ eq '');
		my $data_line = $_ . "\t" . $item . "\n";
		$collective_bed{$data_line} = 1;
	}
	close $temp_stream_processed;
	unlink $outfile_single_phenotype;
	unlink $outfile_single_phenotype_processed;
	#run bedtools, open file, save data in %collective_bed, delete temp file
}
open (my $segment_stream,  q{>}, $outfile_segments) or die("Unable to open $outfile_segments : $!");
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
#open (my $segment_stream,  q{>}, $outfile_segments) or die("Unable to open $outfile_segments : $!");
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