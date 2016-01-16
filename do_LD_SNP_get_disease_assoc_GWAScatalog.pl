#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#THE GWASCATALOG IM USING IS MAPPED TO DBSNP 141 and hg38!! YOU CAN'T DO CHECKS BY CHR IF U USE HAPLOREG (hg19) or snap (hg18) results


#I want to see if any of the candidate snps coming out of do_association_SNPTEST_* are associated in any way with disease or traits in common databases.
#First I input my snps to SNAP/haploreg and get all SNPS in strong LD with my initial snps
#Then I use this script to scan all those SNPS in strong LD with my initial snps with some disease data.

#This script here uses the disease data from the GWAS catalog
#http://www.genome.gov/gwastudies/
#I use a version I postprocess by adding other MS gwas data from the Beecham Publication

#INPUTS:
#-1 SNAP/haploreg output file
#-GWAS catalog
#-2 PVALUE threshold (genome wide significance? -logp = 5?)

#OUTPUT
#-1 list of phenotypes/traits passing the requirements if any

#GWAS CATALOG FIELDs:
#0 Date Added to Catalog
#1 PUBMEDID
#2 First Author
#3 Date
#4 Journal
#5 Link
#6 Study
#7 Disease/Trait
#8 Initial Sample Size
#9 Replication Sample Size
#10 Region
#11 Chr_id
#12 Chr_pos
#13 Reported Gene(s)
#14 Mapped_gene
#15 Upstream_gene_id
#16 Downstream_gene_id
#17 Snp_gene_ids
#18 Upstream_gene_distance
#19 Downstream_gene_distance
#20 Strongest SNP-Risk Allele
#21 SNPs
#22 Merged
#23 Snp_id_current
#24 Context
#25 Intergenic
#26 Risk Allele Frequency
#27 p-Value
#28 Pvalue_mlog
#29 p-Value (text)
#30 OR or beta
#31 95% CI (text)
#32 Platform [SNPs passing QC]
#33 CNV




#structure of snap LD snps to test for disease association
#SNP	Proxy	Distance	RSquared	DPrime	Arrays	Chromosome	Coordinate_HG18
#rs13386439	rs13386439	0	1.000	1.000	None	chr2	212028356
#rs13386439	rs76690544	875	0.877	1.000	None	chr2	212029231
#rs13386439	rs6748422	1195	0.877	1.000	None	chr2	212029551
#rs13386439	rs1992028	1873	0.877	1.000	None	chr2	212026483
#rs13386439	rs2033646	2747	0.877	1.000	AN,A5,A6	chr2	212031103
	
#structure of HAPLOREG LD snps to test
#13 query snps
#chr	pos	r2	D'	is_query_snp	rsID	ref	alt	AFR	AMR	ASN	EUR	GERP_cons	SiPhy_cons	Promoter_ENCODE	Enhancer_ENCODE	Promoter_Roadmap	Enhancer_Roadmap	DNAse	Proteins	eQTL	Motifs	GENCODE_id	GENCODE_name	GENCODE_direction	GENCODE_distance	RefSeq_id	RefSeq_name	RefSeq_direction	RefSeq_distance	dbSNP_functional_annotation
#8	110486627	0.95	0.99	0	rs1457283	A	G	0.35	0.49	0.34	0.34	0	0	.	.	.	ST.MUC,12_EnhWk2	.	.	.	RXRA_disc2	ENSG00000205038.7	PKHD1L1	0	0	NM_177531	PKHD1L1	0	0	INT
									
my $infile_data;
my $source_of_data;
my $PROBABILITY; #make it user selectable?
GetOptions(
		'i=s'        =>\$infile_data,
		'source=s'   =>\$source_of_data,
        'p=f'        =>\$PROBABILITY,
);
my $infile_ph_data = '/net/isi-scratch/giuseppe/indexes/GWAS/gwascatalog_20140218_hg19_PLUS_beechamMS.txt';

#$infile_data = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_diffbind/overlap_3/GAT/all_sorted_plus_LD_haploreg_CEU_0.8.txt";
#$source_of_data = "haploreg";

if(!$infile_data){
	print "USAGE: do_LD_SNP_get_disease_assoc.GWAS.pl -i=<LD_SNPS> -source=[snap|haploreg] -p=<SNP_PVAL>\n";
    print "<LD_SNPS> tab separated output of SNAP/Haploreg for snps proxy search\n";
    print "[snap|haploreg] source of LD info\n";
    print "(optional)<SNP_PVAL> significance level for association (default = 5*10-8)\n";
    exit 1;
}
if(!$source_of_data){
	print "USAGE: do_LD_SNP_get_disease_assoc.GWAS.pl -i=<LD_SNPS> -source=[snap|haploreg] -p=<SNP_PVAL>\n";
    print "<LD_SNPS> tab separated output of SNAP/Haploreg for snps proxy search\n";
    print "[snap|haploreg] source of LD info\n";
    print "(optional)<SNP_PVAL> significance level for association (default = 5*10-8)\n";
    exit 1;
}
#$PROBABILITY = 0.00001 if(!$PROBABILITY);
#Ram says the genome wide threshold is 10-8
#Wikipedia says it's 5x10-8
$PROBABILITY = 0.00000005 if(!$PROBABILITY);

#build a hash from the LD dataset
my %snp_data;
my %snp_coord;
if($source_of_data eq 'snap'){
	#SNP	Proxy	Distance	RSquared	DPrime	Arrays	Chromosome	Coordinate_HG18
	open (my $instream,  q{<}, $infile_data) or die("Unable to open $infile_data : $!");
	while(<$instream>){
		chomp;
		next if($_ =~ /^SNP/);
		next if($_ =~ /WARNING/);
		my ($tag_snp, $rsid, $distance, $rsquared, $chr) = (split /\t/)[0,1,2,3,6];
		#next if($rsquared < 0.9);
		#remove any white spaces from both end of RSid
		#$id =~ s/^\s+|\s+$//g;
		#$snps_names{$rsid} = 1;
		$snp_data{$rsid}{'chr'} = $chr;
		$snp_data{$rsid}{'distance'} = $distance;
		$snp_data{$rsid}{'rsquared'} = $rsquared;
		$snp_data{$rsid}{'tagsnp'} = $tag_snp;
	}
	close $instream;
}elsif($source_of_data eq 'haploreg'){
	#first find the tag snp (ie the set of those YOU provided haploreg with)
	#chr	pos	r2	D'	is_query_snp	rsID	ref	alt
	open (my $instream,  q{<}, $infile_data) or die("Unable to open $infile_data : $!");
	while(<$instream>){
		next if($_ =~ /query snp/);
		next if($_ =~ /^chr/);
		next if($_ eq '');
		my ($chr, $pos, $rsquared, $is_tag_snp, $rsid) = (split /\t/)[0,1,2,4,5];
		#$snps_names{$rsid} = 1;
		$snp_data{$rsid}{'chr'} = $chr;
		$snp_data{$rsid}{'distance'} = $pos; #this is actually a position!!
		$snp_data{$rsid}{'rsquared'} = $rsquared;
		$snp_data{$rsid}{'tagsnp'} = $is_tag_snp;
		my $coord = $chr . '-' . $pos;
		$snp_coord{$coord} = 1;
	}
	close $instream;
}else{
	print "ERROR: -source field: $source_of_data not recognised. Aborting..\n";
	exit -1;
}


#open and process the GWAS catalog data
open (my $instream,  q{<}, $infile_ph_data) or die("Unable to open $infile_ph_data : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^Date/);
	#check the rs_id. If it's in the hash then check the pvalue. If it's below threshold then print a summary line on stdout?
	
	my ($phenotype, $chr, $pos, $rsID, $pvalue) =  (split /\t/)[7, 11, 12, 21,27];
	next if(!$rsID);
	next if($rsID =~ /NR/);
	#what about multiple rsids?
	#eg rs204999,rs9268528,rs9268542,rs6903608,rs2858870
	#split and test all of them
	my @rsID = split("," , $rsID);
	
	foreach my $ID (@rsID){
		if($snp_data{$rsID}){
			if($pvalue  <= $PROBABILITY){
				print	STDOUT 
					$snp_data{$rsID}{'chr'} . "\t" .
					$snp_data{$rsID}{'distance'} . "\t" .
				 	$rsID . "\t" .
					$snp_data{$rsID}{'tagsnp'}   . "\t" .
					$snp_data{$rsID}{'rsquared'} . "\t" .
					$_ . "\n";
			}
		}			
	}
}
close $instream;


#uncomment this if you are sure you use all hg19 or all hg38 

#if the data is from Haploreg, test via coordinates (rsID change)
if($source_of_data eq 'haploreg'){
	print "\nTESTING BY COORDINATES..\n";
	open (my $instream,  q{<}, $infile_ph_data) or die("Unable to open $infile_ph_data : $!");
	while(<$instream>){
		chomp;
		next if($_ =~ /^Date/);
		my ($phenotype, $chr, $pos, $rsID, $pvalue) =  (split /\t/)[7, 11, 12, 21,27];
		next if(!$chr);
		next if ($chr eq '');

		my $coord =  $chr . '-' . $pos;
		
		if($snp_coord{$coord} && ($pvalue  <= $PROBABILITY)  ){
			print	STDOUT 
					$snp_data{$rsID}{'chr'} . "\t" .
					$snp_data{$rsID}{'distance'} . "\t" .
				 	$rsID . "\t" .
					$snp_data{$rsID}{'tagsnp'}   . "\t" .
					$snp_data{$rsID}{'rsquared'} . "\t" .
					$_ . "\n";
		}
	}
	close $instream;
}
print "FINISHED\n";