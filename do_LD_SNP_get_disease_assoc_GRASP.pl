#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#I want to see if any of the candidate snps coming out of do_association_SNPTEST_* are associated in any way with disease or traits in common databases.
#First I input my snps to SNAP or HAPLOREG (seems to have more SNP ids) and get all SNPS in LD with my initial snps
#Then I use this script to scan all those SNPS in  LD with my initial snps with some disease data.

#This script here uses the disease data from GRASP
#paper: 
#GRASP: analysis of genotype-phenotype results from 1390 genome-wide association studies and corresponding open access database
#http://bioinformatics.oxfordjournals.org/content/30/12/i185.full 

#data:
#http://apps.nhlbi.nih.gov/Grasp/Updates.aspx

#INPUTS:
#-1 SNAP output file OR HAPLOREG output file
#-2 which of the two it is
#-3 PVALUE threshold (genome wide significance? -logp = 5?)

#OUTPUT
#-1 list of phenotypes/ traits passing the requirements if any

#structure of SNAP LD snps to test for disease association
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


#structure of disease data
#NHLBIkey	HUPfield	LastCurationDate	CreationDate	SNPid(dbSNP134)	chr(hg19)	pos(hg19)	PMID	SNPid(in paper)	LocationWithinPaper	Pvalue	Phenotype	PaperPhenotypeDescription	PaperPhenotypeCategories	DatePub	InNHGRIcat(as of 3/31/12)	Journal	Title	IncludesMale/Female Only Analyses	Exclusively Male/Female	Initial Sample Description	Replication Sample Description	Platform [SNPs passing QC]	GWASancestryDescription	TotalSamples(discovery+replication)	TotalDiscoverySamples	European Discovery	African Discovery	East Asian Discovery	Indian/South Asian Discovery	Hispanic Discovery	Native Discovery	Micronesian Discovery	Arab/ME Discovery	Mixed Discovery	Unspecified Discovery	Filipino Discovery	Indonesian Discovery	Total replication samples	European Replication	African Replication	East Asian Replication	Indian/South Asian Replication	Hispanic Replication	Native Replication	Micronesian Replication	Arab/ME Replication	Mixed Replication	Unspecified Replication	Filipino Replication	Indonesian Replication	InGene	NearestGene	InLincRNA	InMiRNA	InMiRNABS	dbSNPfxn	dbSNPMAF	dbSNPalleles/het/se	dbSNPvalidation	dbSNPClinStatus	ORegAnno	ConservPredTFBS	HumanEnhancer	RNAedit	PolyPhen2	SIFT	LS-SNP	UniProt	EqtlMethMetabStudy
#203831461	Jan2014	8/17/12	8/17/12	3	13	32446842	20383146	rs3	Full Scan	0.034	Chronic kidney disease	Chronic kidney disease (CKD) and renal traits	Renal;Chronic kidney disease;Quantitative trait(s)	4/11/2010	y	Nat Genet	New loci associated with kidney function and chronic kidney disease.	n	n	Up to 67093 EA individuals	Up to 22982 EA individuals	Affymetrix & Illumina [~2.5 million] (imputed)	European	90075	67093	67093												22982	22982												(EEF1DP3)					Intron	T;0.082	C/T;0.134821;0.221887	YES										
#204538422	Jan2014	8/17/12	8/17/12	3	13	32446842	20453842	rs3	FullScan	0.0209937	Rheumatoid arthritis	Rheumatoid arthritis	Inflammation;Arthritis;Rheumatoid arthritis	5/9/2010	y	Nat Genet	Genome-wide association study meta-analysis identifies seven new rheumatoid arthritis risk loci.	n	n	5539 EA cases, 20169 EA controls	6768 EA cases, 8806 EA controls	Affymetrix & Illumina [~2716259] (imputed)	European	41282	25708	25708												15574	15574												(EEF1DP3)					Intron	T;0.082	C/T;0.134821;0.221887	YES										

my $infile_data;
my $source_of_data;
my $PROBABILITY; #make it user selectable?
GetOptions(
		'i=s'        =>\$infile_data,
		'source=s'   =>\$source_of_data,
        'p=i'        =>\$PROBABILITY,
);

if(!$infile_data){
	print "USAGE: do_LD_SNP_get_disease_assoc.GRASP.pl -i=<LD_SNPS> -source=[snap|haploreg] -p=<SNP_PVAL>\n";
    print "<LD_SNPS> tab separated output of SNAP/Haploreg for snps proxy search\n";
    print "[snap|haploreg] source of LD info\n";
    print "(optional)<SNP_PVAL> significance level for association (default = 5*10-8)\n";
    exit 1;
}
if(!$source_of_data){
	print "USAGE: do_LD_SNP_get_disease_assoc.GRASP.pl -i=<LD_SNPS> -source=[snap|haploreg] -p=<SNP_PVAL>\n";
    print "<LD_SNPS> tab separated output of SNAP/Haploreg for snps proxy search\n";
    print "[snap|haploreg] source of LD info\n";
    print "(optional)<SNP_PVAL> significance level for association (default = 5*10-8)\n";
    exit 1;
}

my $infile_ph_data = '/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/Jan9graspForNCBI';


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
	#getting the query snp is a mess. you have cases like
	#5	139340779	1	1	1	rs6580323
	#5	38851925	0.97	0.99	0	rs357289
	#5	38855345	1	1	1	rs548287
	#5	38857366	0.93	0.99	0	rs420444
	#so same chromosome, two query snps. You need to check the position somehow
	
	my %query_snps;
	my $counter = 1;
	#first pass: get the list of query snps ORDERED by counter in the way they appear in the data
	open (my $instream,  q{<}, $infile_data) or die("Unable to open $infile_data : $!");
	while(<$instream>){
		next if($_ =~ /query snps/);
		next if($_ =~ /^chr/);
		next if($_ eq '');
		my ($is_query_snp, $rsid) = (split /\t/)[4,5];
		if($is_query_snp == 1){
			$query_snps{$counter} = $rsid;
			$counter++;
		}
	}
	close $instream;	
	
	$counter = 0; #this needs to be incremented at the transition between ld groups
	my $chr_pos_buffer = '';
	#second pass
	#get data for all snps but label them with their query snp
	#chr	pos	r2	D'	is_query_snp	rsID	ref	alt
	open ($instream,  q{<}, $infile_data) or die("Unable to open $infile_data : $!");
	while(<$instream>){
		next if($_ =~ /query snps/);
		next if($_ =~ /^chr/);
		my ($chr, $pos, $rsquared, $is_tag_snp, $rsid) = (split /\t/)[0,1,2,4,5];
		my $pos_first = substr($pos, 0, 1);
		my $chr_pos_sig = $chr . '_' . $pos_first;
		
		$counter++ if($chr_pos_sig ne $chr_pos_buffer);
		$chr_pos_buffer = $chr_pos_sig;

		$snp_data{$rsid}{'chr'} = $chr;
		$snp_data{$rsid}{'distance'} = $pos; #this is actually a position!!
		$snp_data{$rsid}{'rsquared'} = $rsquared;
		$snp_data{$rsid}{'tagsnp'} = $query_snps{$counter};
		my $hg19_coord = $chr . '-' . $pos;
		$snp_coord{$hg19_coord} = 1;	
	}
	close $instream;
}else{
	print "ERROR: -source field: $source_of_data not recognised. Aborting..\n";
	exit -1;
}


#disease catalogue
#build header
#compose data from the LD haploreg stuff and disease stuff
#this comes from the GRASP catalogue
#"NHLBIkey\tHUPfield\tLastCurationDate\tCreationDate\tSNPid(dbSNP134)\tchr(hg19)\tpos(hg19)\tPMID\tSNPid(in paper)\tLocationWithinPaper\tPvalue\tPhenotype\tPaperPhenotypeDescription\tPaperPhenotypeCategories\tDatePub\tInNHGRIcat(as of 3/31/12)\tJournal\tTitle\tIncludesMale/Female Only Analyses\tExclusively Male/Female\tInitial Sample Description\tReplication Sample Description\tPlatform [SNPs passing QC]\tGWASancestryDescription\tTotalSamples(discovery+replication)\tTotalDiscoverySamples\tEuropean Discovery\tAfrican Discovery\tEast Asian Discovery\tIndian/South Asian Discovery\tHispanic Discovery\tNative Discovery\tMicronesian Discovery\tArab/ME Discovery\tMixed Discovery\tUnspecified Discovery\tFilipino Discovery\tIndonesian Discovery\tTotal replication samples\tEuropean Replication\tAfrican Replication\tEast Asian Replication\tIndian/South Asian Replication\tHispanic Replication\tNative Replication\tMicronesian Replication\tArab/ME Replication\tMixed Replication\tUnspecified Replication\tFilipino Replication\tIndonesian Replication\tInGene\tNearestGene\tInLincRNA\tInMiRNA\tInMiRNABS\tdbSNPfxn\tdbSNPMAF\tdbSNPalleles/het/se\tdbSNPvalidation\tdbSNPClinStatus\tORegAnno\tConservPredTFBS\tHumanEnhancer\tRNAedit\tPolyPhen2\tSIFT\tLS-SNP\tUniProt\tEqtlMethMetabStudy";


my $HEADER = "CHR\tPOS\trsID\tQUERY_SNP\tRSQUARED\tNHLBIkey\tHUPfield\tLastCurationDate\tCreationDate\tSNPid(dbSNP134)\tchr(hg19)\tpos(hg19)\tPMID\tSNPid(in paper)\tLocationWithinPaper\tPvalue\tPhenotype\tPaperPhenotypeDescription\tPaperPhenotypeCategories\tDatePub\tInNHGRIcat(as of 3/31/12)\tJournal\tTitle\tIncludesMale/Female Only Analyses\tExclusively Male/Female\tInitial Sample Description\tReplication Sample Description\tPlatform [SNPs passing QC]\tGWASancestryDescription\tTotalSamples(discovery+replication)\tTotalDiscoverySamples\tEuropean Discovery\tAfrican Discovery\tEast Asian Discovery\tIndian/South Asian Discovery\tHispanic Discovery\tNative Discovery\tMicronesian Discovery\tArab/ME Discovery\tMixed Discovery\tUnspecified Discovery\tFilipino Discovery\tIndonesian Discovery\tTotal replication samples\tEuropean Replication\tAfrican Replication\tEast Asian Replication\tIndian/South Asian Replication\tHispanic Replication\tNative Replication\tMicronesian Replication\tArab/ME Replication\tMixed Replication\tUnspecified Replication\tFilipino Replication\tIndonesian Replication\tInGene\tNearestGene\tInLincRNA\tInMiRNA\tInMiRNABS\tdbSNPfxn\tdbSNPMAF\tdbSNPalleles/het/se\tdbSNPvalidation\tdbSNPClinStatus\tORegAnno\tConservPredTFBS\tHumanEnhancer\tRNAedit\tPolyPhen2\tSIFT\tLS-SNP\tUniProt\tEqtlMethMetabStudy";
open (my $instream,  q{<}, $infile_ph_data) or die("Unable to open $infile_ph_data : $!");
print STDOUT $HEADER, "\n";
while(<$instream>){
	chomp;
	next if($_ =~ /^NHLBIkey/);
	#check the rs_id. If it's in the hash then check the pvalue. If it's below threshold then print a summary line on stdout?
	#rsId: field 4 (hg19)
	#chr: field 5
	#pos:field 6
	#pvalue: field 10
	#phenotype: field 11
	
	my ($rsID, $chr, $pos, $pvalue, $phenotype, $ph_description, $ph_category) = (split /\t/)[4,5,6,10,11,12,13];
	next if(!$pvalue);
	
	#REMOVE THIS?
	next if($phenotype =~ /^Gene expression of/); #remove qtl
	#REMOVE THIS?
	
	$rsID = 'rs' . $rsID;
	if($snp_data{$rsID}){
		if($pvalue  <= $PROBABILITY){
			print	STDOUT 
					$snp_data{$rsID}{'chr'} . "\t" .
					$snp_data{$rsID}{'distance'} . "\t" .
				 	$rsID . "\t" .
					$snp_data{$rsID}{'tagsnp'}   . "\t" .
					$snp_data{$rsID}{'rsquared'} . "\t" .
					$_ . "\n";

			#print STDOUT $_, "\n\n";
		}
	}
}
close $instream;

#if the data is from Haploreg, test via coordinates (rsID change)
if($source_of_data eq 'haploreg'){
	print STDOUT "\nTESTING BY COORDINATES..\n";
	open (my $instream,  q{<}, $infile_ph_data) or die("Unable to open $infile_ph_data : $!");
	while(<$instream>){
		chomp;
		next if($_ =~ /^NHLBIkey/);
		my ($rsID, $chr, $pos, $pvalue, $phenotype, $ph_description, $ph_category) = (split /\t/)[4,5,6,10,11,12,13];
		
		#REMOVE THIS?
		next if($phenotype =~ /^Gene expression of/); #remove qtl
		#next if($phenotype =~ /^Methylation levels /); #remove qtl
		#REMOVE THIS?
		
		next if(!$pvalue);
		next if(!$chr);
		next if ($chr eq '');
		my $coord = $chr . '-' . $pos;
		$rsID = 'rs' . $rsID;
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