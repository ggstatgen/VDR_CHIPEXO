#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#get all the rsID for the positions in the "interestinghet_consensus.txt" output of do_alleleseq_pp_01_intersect_results.pl
#get the genotypes for that position from all the samples 
#THIS VERSION OUTPUTS VCF files AND USES A HUGE HASH TO STORE THE ORIGINAL VCF DATA. ONLY RUN THIS ON HIGH MEMORY MACHINES.

#usage scenarios:
#you might need the output of this function to send them to SNAP/Haploreg, get the snps in LD, check for disease association with the do_SNAP_... scripts, and filter via bedtools to only get snps in motifs.

#There is the problem that 15 samples have 1kg data and 5 only imputed from 1kg.
#the best way I have found to address this is the following:
#1 filter genotypes for the 15 samples from the 1kg vcf using do_vcf_slice - get vcf.gz
#2 filter genotypes for the 5 samples from the IMPUTE2 vcf using do_vcf_slice - get vcf.gz
#3 index with tabix
#4 merge the two resulting vcfs in a comprehensive vcf.gz with genotypes for the 20 samples using /net/isi-scratch/giuseppe/tools/HTSLIB/bcftools
#use this for the operations in this file
my $INPUT_VARIANTS = '/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/ALLELESEQ_20samples_merged.vcf.gz';
#my $INPUT_VARIANTS = '/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/test.vcf.gz';

#new vcf headers, modify if needed
my $INFO_PK = 'ASB_PK';
my $INFO_SA = 'ASB_SA';
my $header_alleleseq_pk = "##INFO=<ID=$INFO_PK,Number=0,Type=Flag,Description=\"event in VDR CPo3 peak (alleleseq)\">";
my $header_alleleseq_sa = "##INFO=<ID=$INFO_SA,Number=1,Type=String,Description=\"samples testing positive for allelic imbalance (alleleseq)\">";

#structure:
###fileformat=VCFv4.1
###FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype calls">
###FORMAT=<ID=GP,Number=3,Type=Float,Description="Genotype call probabilities">
##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA06986 NA06989 NA06997 NA07029 NA07045 NA10831 NA10846 NA10847 NA11829 NA11832 NA11918 NA11919 NA12264 NA12383 NA12489 NA12716 NA12752 NA12872 NA19189 NA19190 NA19191
#chr1    10583   rs58108140      G       A       .       .       .       GT:GP   0/1:0.031,0.963,0.006   0/1:0.032,0.961,0.006   ./.:0.733,0.226,0.041   ./.:0.784,0.206,0.01    ./.:0.468,0.393,0.139   ./.:0.755,0.24,0.005    ./.:0.841,0.1

#OUTPUT
#vcf with entries
my $infile_alleleseq;
my $allsites; #flag. If set, gets ASB sites even when they don't fall in a peak bed interval
my $ucsc_format; #flag. If set, use 'hg19' chromosome names instead of b37 in output (useful to post-process variants using Funseq)
GetOptions(
        'i=s'        =>\$infile_alleleseq,
        'ucsc_chr'   =>\$ucsc_format,
        'allsites'   =>\$allsites
);

#$infile_alleleseq = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/CPO3_PEAKS/OVERLAP/interestingHets_consensus_inpeaks_o2.txt";

if(!$infile_alleleseq){
	print "USAGE: do_alleleseq_pp_03_ASB_to_VCF.pl -i=<INFILE_ASSOC> -ucsc_chr -allsites\n";
	print "<INFILE_ASSOC> file interestingHets_overlaps.txt from do_alleleseq_pp_02..\n";
	print "(opt)<ucsc_chr> flag; if set, save chromosome names in UCSC format (chrX). Default: b37 format (X)\n";	
	print "(opt)<allsites> flag; if set, retrieve asb sites even when they did not fall in peak interval file.\n";
    exit 1;
}
my $TOOL_PATH       = '/net/isi-scratch/giuseppe/tools';
my $BCFTOOLS        = $TOOL_PATH . '/HTSLIB/bcftools/bcftools';
my $VCF_SORT        = $TOOL_PATH . '/vcftools_0.1.12b/bin/vcf-sort';
#my $BGZIP           = $TOOL_PATH . 'HTSLIB/htslib/bgzip';
#my $TABIX           = $TOOL_PATH . '/HTSLIB/htslib/tabix';

my ($filename, $directory) = fileparse($infile_alleleseq);
$filename =~ s/(.*)\..*/$1/;

my $out_file; my $out_file_sorted;
if($ucsc_format){
	$out_file = $directory . $filename . '_hg19.vcf';
	$out_file_sorted = $directory . $filename . '_sorted_hg19.vcf';	
}else{
	$out_file = $directory . $filename . '.vcf';
	$out_file_sorted = $directory . $filename . '_sorted.vcf';
}
#my $out_file_sorted_gzipped = $directory . $filename . '_sorted.vcf.gz';

#process header------------------------------------
my $header = `$BCFTOOLS view -h $INPUT_VARIANTS`;;
#do some clean up? remove all the "contigs" entries
my @header = split("\n", $header);
my @out_header;
foreach my $item (@header){
	push @out_header, $item unless($item =~ /\#\#contig=(.*)/ || $item =~ /\#\#bcftools_(.*)/);
}
#get first and last row
my $header_vcf_spec = shift @out_header;
my $header_column_names = pop @out_header;
#add new rows
push @out_header, $header_alleleseq_pk;
push @out_header, $header_alleleseq_sa;
#sort
my @sorted_out_header = sort @out_header;
#re-add top and bottom row
unshift @sorted_out_header, $header_vcf_spec;
push @sorted_out_header, $header_column_names;
my $out_header = join("\n", @sorted_out_header);
#process header-----------------------------------

#############
#HASH containing variants
#############
my %hash_variants;
my %hash_variants_count;
tie *FILE,   'IO::Zlib', $INPUT_VARIANTS, "rb";
while (<FILE>)	{ 
	chomp;
	next if($_ eq '');
	next if($_ =~ /^##/);
	
	my @fields = split("\t", $_);
	
	#the key will be "chr-position"
	my $key = $fields[0] . '-' . $fields[1];
	
	$hash_variants{$key} = $_;
	if(!$hash_variants_count{$key}){
		$hash_variants_count{$key} = 1;
	}else{
		$hash_variants_count{$key} += 1;
	}
}
close FILE;
#This should have taken a long time. Once it's done, you can map the asb events to rsID and build an output VCF
#map the actual instances, remove associations not in binding sites (unless allsites), output enriched VCF files, then merge them


my %vcf_data;
open (my $instream,  q{<}, $infile_alleleseq) or die("Unable to open $infile_alleleseq : $!");
while(<$instream>){
	my $vcf_line;
	chomp;
	
	#checks-----------
	next if($_ eq '');
	next if($_ =~ /^chrm/);
	my ($a_chr,$a_snppos,$a_bindingSite) = (split /\t/)[0,1,3];
	next if(!$a_chr);
	next if(!$a_snppos);
	#binding sites yes/no?
	unless($allsites){
		next unless ($a_bindingSite =~ /1/);
	}
	#checks-----------
	my @vcf_line;
	#now get all fields
	my $alleleseq_line = $_;
	my @fields = split("\t", $_);
	#test chr and location against vcf hash

	my $key = $a_chr . '-' . $a_snppos;
	if($hash_variants{$key}){
		if($hash_variants_count{$key} > 1){
			#IF we are here, there are two or more variants at the same position.
			#in my experience, this is when an indel overlaps a snp, where the indel is on one allele and the snp on the other
			#When vcf2diploid builds the personalised genomes, it SKIPS these. Therefore these variants are not incorporated 
			#in the personalised genomes
			print STDERR "key $key is associated to more than one vcf entry. Skipping\n"; # you could check the reference and the alternate?
			next;
		}
		#subroutine to take vcf data and alleleseq data and build a new combined vcf line
		my $vcf_line  = get_new_vcf_line($alleleseq_line, $hash_variants{$key});
		#save these lines in a hash, then write file when you're done. Or write immediately?
		$vcf_data{$vcf_line} = 1;
	}else{
		print STDERR "ASB event: $alleleseq_line - no matching rsID found for this position in vcf reference data. Skipping..\n";
		next;
	}
}
close $instream;

#write an output vcf, including a header describing the new fields
open (my $outstream,  q{>}, $out_file) or die("Unable to open $out_file : $!");
print $outstream $out_header, "\n"; #taken from the reference data
foreach my $item (keys %vcf_data){ print $outstream $item . "\n"; }
close $outstream;
#sort with vcf-tools, compress with bgzip, index with tabix
system "$VCF_SORT -c $out_file > $out_file_sorted";
#system "$BGZIP $out_file_sorted";
#system "$TABIX -p vcf $out_file_sorted_gzipped";
unlink $out_file;
system "mv $out_file_sorted $out_file";

############
#SUBS
############

#sample alleleseq input line:
#14      107169846       C       1       8       NA06986,NA06989,NA10847,NA11829,NA12489,NA19189,NA19213,NA19248
#0 - chromosome
#1 - position (hg19)
#2 - reference
#3 - in cp_o3 peak (1/0)
#4 - number of samples where asb event detected
#5 - sample list where asb detected

#sample VCF 4.1
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT [GENOTYPES]
#1       10583   rs58108140      G       A       100     PASS    AA=.;AC=6;AF=0.14;AFR_AF=0.04;AMR_AF=0.17;AN=32;ASN_AF=0.13;AVGPOST=0.7707;ERATE=0.0161;EUR_AF=0.21;LDAF=0.2327;RSQ=0.4319;SNPSOURCE=LOWCOV;THETA=0.0046;VT=SNP GT:DS:GL
#0 - chromosome
#1 - position (hg19)
#2 - rsid
#3 - ref
#4 - alt
#5 - qual
#6 - filter
#7 - info
#8 - format
#9..x - genotypes
sub get_new_vcf_line{
	my ($alleleseq_line, $vcf_line) = @_;
	my $info_field;
	#create info line from alleleseq fields number 3,5
	#split vcf data line, get base info and genotypes
	
	my @vcf_line = split("\t", $vcf_line);	
	my @alleleseq_line = split("\t", $alleleseq_line);
	
	if($alleleseq_line[3] eq '1'){
		$info_field = $vcf_line[7] . ';' . $INFO_PK . ';' . $INFO_SA . '=' . $alleleseq_line[5];
	}else{
		$info_field = $vcf_line[7] . ';' . $INFO_SA . '=' . $alleleseq_line[5];
	}
	#remove the old info field from @vcf_line and plug in the new one
	splice @vcf_line, 7, 1; #remove old info field
	splice @vcf_line, 7, 0, $info_field; #insert new info field
	
	#if UCSC output is desired, replace $vcf_line[0] with 'chr' . $vcf_line[0]
	$vcf_line[0] = 'chr' . $vcf_line[0] if($ucsc_format);
	
	return join("\t", @vcf_line);
}
