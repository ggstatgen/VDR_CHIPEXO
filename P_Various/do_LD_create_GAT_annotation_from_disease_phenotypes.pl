#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use IO::Zlib;

#once you've found that some VDR-BV/hcBV/ucBV snps are in strong LD with GWAS snps for some trait/disease, you want to assess the significance of this association

#This file accepts as input an LD_GWAS interval file:
#/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/LD_PLINK_CEU_GRASP2final.vcf.gz
#..
#/net/isi-scratch/giuseppe/indexes/GWAS_GRASP/LD_PLINK_CEU_GRASP2_plus_Beecham2013.bed.gz

#and a txt file with a list of phenotypes (a subset of those present in the file above) including only those you wish to test significance against.
#this will be similar to  do_GRASPcatalog_get_bed_for_GAT_byphenotype.pl, but there I was getting the GWAS snps for the designated phenotypes (not the LD intervals containing them)

#output: a bed file (one per catalog), each row labelled by phenotype

my $phenotype_file;
my $catalog_file;
my $LUMP_PHENOTYPES;  #if lump phenotypes, I will do a reg ex test on the phenotype match. So things like
#Rheumatoid arthritis
#Rheumatoid arthritis (ACPA-positive)
#Rheumatoid arthritis and celiac disease
#Rheumatoid arthritis, combined control dataset
#Rheumatoid arthritis, cyclic citrullinated peptide (CCP) positive

#should all go as "Rheumatoid arthritis"
#TODO

GetOptions(
        'i=s'       =>\$catalog_file,
        'ph=s'      =>\$phenotype_file #,
        #'lump'      =>\$LUMP_PHENOTYPES
);
if(!$catalog_file){
     print "USAGE: do_LD_create_GAT_annotation_from_disease_phenotype_intervals.pl -i=<LD_CAT> -ph=<PHENOTYPES>\n";
     print "<LD_CAT> catalog of PLINK LD blocks overlapping GRASP/GWAS SNPs\n";
     print "<PHENOTYPES> text file with one phenotype per row, including only the phenotypes you want included in the output .bed (phenotypes should be compatibe with LD_CAT)\n";
     print "NEEDS IMPROVEMENT! IT WILL ONLY FIND EXACT PHENOTYPE NAMES\n";
     exit 1;
}
if(!$phenotype_file){
     print "USAGE: do_LD_create_GAT_annotation_from_disease_phenotype_intervals.pl -i=<LD_CAT> -ph=<PHENOTYPES>\n";
     print "<LD_CAT> catalog of PLINK LD blocks overlapping GRASP/GWAS SNPs with p < 5 * 10^8\n";
     print "<PHENOTYPES> text file with one phenotype per row, including only the phenotypes you want included in the output .bed (phenotypes should be compatibe with LD_CAT)\n";
     print "NEEDS IMPROVEMENT! IT WILL ONLY FIND EXACT PHENOTYPE NAMES\n";
     exit 1;
}
my $TOOL_PATH       = '/net/isi-scratch/giuseppe/tools';
my $BEDTOOLS  = $TOOL_PATH . '/bedtools-2.22.1/bin';

#I want the output to be in the directory of the phenotype file
my($basename, $directory) = fileparse($phenotype_file);
$basename =~ s/(.*)\..*/$1/;
my $outfile = $directory . "\/" .  $basename . '_segments.bed';

#get list
my %phenotype_map;
open (my $instream,  q{<}, $phenotype_file) or die("Unable to open $phenotype_file : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	$phenotype_map{lc($_)} = 1;
}
close $instream;


#the catalog, structure-wise, is a [bed_ldblock] + [vcf_gwas_snp] TSV file
#chr	ldstart	ldstop	snpsinld	kb_block	chr_gwas_snp	pos_gwas_snp	id_gwas_snp	-	-	-	phenotype	phen_desc	1
#you need 0,1,2,11 (ld info (0,1,2) + associated phenotype (11))
#FOR EACH phenotype, you need to merge all overlapping LD blocks. Whilst the same block might be listed multiple times, for different phenotypes, you want to see it once per phenotype. Therefore, create a temp bed, one per phenotype, merge overlapping blocks, join into a final bed


#create 1bed for each phenotype
#merge overlapping intervals in each bed
#create final bed
#create one bed for each of the terms above, if the term has > $min_snp_number
my %collective_bed;
foreach my $item (sort keys %phenotype_map){
	my $outfile_single_phenotype           = $directory . "\/" .  $basename . '_phenotype_' . '_temp.bed';
	my $outfile_single_phenotype_processed = $directory . "\/" .  $basename . '_phenotype_' . '_temp_out.bed';
	tie *FILE,   'IO::Zlib', $catalog_file, "rb";
	open (my $temp_stream,  q{>}, $outfile_single_phenotype) or die("Unable to open $outfile_single_phenotype : $!");
	while(<FILE>){
		chomp;
		next if($_ eq '');
		my ($ld_block_chr, $ld_block_start, $ld_block_stop, $phenotype) = (split /\t/)[0,1,2,11];
		if(lc($phenotype) eq $item){
			my $line = $ld_block_chr . "\t" . $ld_block_start . "\t" . $ld_block_stop  . "\t" . $item;
			print $temp_stream 	$line, "\n";
		}
	}
	close FILE;
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
}
open (my $segment_stream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
foreach my $item (sort keys %collective_bed){ print $segment_stream $item; }
close $segment_stream;