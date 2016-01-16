#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#I want to see if any of the candidate snps coming out of do_association_SNPTEST_* are associated in any way with disease or traits in common databases.
#First I input my snps to SNAP and get all SNPS in strong LD with my initial snps
#Then I use this script to scan all those SNPS in strong LD with my initial snps with some disease/trait data.

#This script here uses the trait data from the Pritchard eQTL collection, LCL only
#http://eqtl.uchicago.edu/Home.html
#the file is a gff v.3 as follows (hg18)

#chr1	Degner2012_dsQTL	Degner_dsQTL	801099	801099	3.68973163336755	.	.	DegnerDSQTL "chr1.801099 a dsQTL for 802000 to 802100" ; Note "Acts in cis "; Plot=<a href="http://eqtl.uchicago.edu/dsQTL_data/FIGURES/caQTL_7_FEB_2012/DegnerDSQTL.chr1.801099.html" >See QTL plots and view/leave comments</a>
#chr1	Degner2012_dsQTL	Degner_dsQTL	846446	846446	5.69615622511135	.	.	DegnerDSQTL "chr1.846446 a dsQTL for 846400 to 846500" ; Note "Acts in cis "; Plot=<a href="http://eqtl.uchicago.edu/dsQTL_data/FIGURES/caQTL_7_FEB_2012/DegnerDSQTL.chr1.846446.html" >See QTL plots and view/leave comments</a>
#chr1	Degner2012_dsQTL	Degner_dsQTL	901912	901912	3.94005811193805	.	.	DegnerDSQTL "chr1.901912 a dsQTL for 901300 to 901400" ; Note "Acts in cis "; Plot=<a href="http://eqtl.uchicago.edu/dsQTL_data/FIGURES/caQTL_7_FEB_2012/DegnerDSQTL.chr1.901912.html" >See QTL plots and view/leave comments</a>
#chr1	Degner2012_dsQTL	Degner_dsQTL	901458	901458	9.22112552799726	.	.	DegnerDSQTL "chr1.901458 a dsQTL for 901400 to 901500" ; Note "Acts in cis "; Plot=<a href="http://eqtl.uchicago.edu/dsQTL_data/FIGURES/caQTL_7_FEB_2012/DegnerDSQTL.chr1.901458.html" >See QTL plots and view/leave comments</a>
#chr1	Degner2012_dsQTL	Degner_dsQTL	904803	904803	4.05340037498487	.	.	DegnerDSQTL "chr1.904803 a dsQTL for 905000 to 905100" ; Note "Acts in cis "; Plot=<a href="http://eqtl.uchicago.edu/dsQTL_data/FIGURES/caQTL_7_FEB_2012/DegnerDSQTL.chr1.904803.html" >See QTL plots and view/leave comments</a>


#INPUTS:
#-1 SNAP output file
#-eqtl pritchard data

#OUTPUT
#-1 list of phenotypes/traits passing the requirements if any


#structure of snap LD snps to test for disease association
#SNP	    Proxy	    Distance	RSquared	DPrime	Arrays	Chromosome	Coordinate_HG18
#rs13386439	rs13386439	0	1.000	1.000	None	chr2	212028356
#rs13386439	rs76690544	875	0.877	1.000	None	chr2	212029231
#rs13386439	rs6748422	1195	0.877	1.000	None	chr2	212029551
#rs13386439	rs1992028	1873	0.877	1.000	None	chr2	212026483
#rs13386439	rs2033646	2747	0.877	1.000	AN,A5,A6	chr2	212031103
									

my $infile_test_data;
my $infile_ph_data;
GetOptions(
		'i=s'        =>\$infile_test_data,

);
$infile_ph_data = '/net/isi-scratch/giuseppe/indexes/GWAS/PRITCHARD_EQTL/LCL.individual.tracks.v3.hg18.gff';

if(!$infile_test_data){
	print "USAGE: do_SNAP_get_disease_assoc_snps_PRITCHARDdnaseqtl.pl -i=<SNAP_SNPS>\n";
    print "<SNAP_SNPS> tab separated output of SNAP for snps proxy search\n";
    print "ATTENTION: Script uses coordinates - it assumes both SNAP and PritchardQTL use hg18\n";
    exit 1;
}


#build a hash from the SNAP dataset
#SNP	Proxy	Distance	RSquared	DPrime	Arrays	Chromosome	Coordinate_HG18
my %SNAP_tagsnp_to_proxysnp;
my %proxy_snps_names;
open (my $instream,  q{<}, $infile_test_data) or die("Unable to open $infile_test_data : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^SNP/);
	next if($_ =~ /WARNING/);
	my ($tag_snp, $proxy_snp, $distance, $rsquared, $chr, $coord_hg18) = (split /\t/)[0,1,2,3,6,7];
	#next if($rsquared < 0.9);
	#remove any white spaces from both end of RSid
	#$id =~ s/^\s+|\s+$//g;
	$proxy_snps_names{$proxy_snp} = 1;
	$SNAP_tagsnp_to_proxysnp{$proxy_snp}{'distance'} = $distance;
	$SNAP_tagsnp_to_proxysnp{$proxy_snp}{'rsquared'} = $rsquared;
	$SNAP_tagsnp_to_proxysnp{$proxy_snp}{'chr'} = $chr;
	$SNAP_tagsnp_to_proxysnp{$proxy_snp}{'tagsnp'} = $tag_snp;
	$SNAP_tagsnp_to_proxysnp{$proxy_snp}{'coord_hg18'} = $coord_hg18;
}
close $instream;

#work on phenotype data
open ($instream,  q{<}, $infile_ph_data) or die("Unable to open $infile_ph_data : $!");
while(<$instream>){
	chomp;
	my ($chr, $start, $stop, $score, $ph_description) = (split /\t/)[0,3,4,5,8];
	#try to overlap coordinates (both hg18)
	#create small test interval
	#TODO TEMP - GET THIS DONE PROPERLY
	$start -= 2;
	$stop  += 2;
	
	foreach my $proxy_id (sort keys %SNAP_tagsnp_to_proxysnp){
		if($SNAP_tagsnp_to_proxysnp{$proxy_id}{'chr'} eq $chr){
			if( ($SNAP_tagsnp_to_proxysnp{$proxy_id}{'coord_hg18'} > $start)  &&  ($SNAP_tagsnp_to_proxysnp{$proxy_id}{'coord_hg18'} < $stop)  ){
				print $proxy_id . "\t" 
								. $SNAP_tagsnp_to_proxysnp{$proxy_id}{'tagsnp'}   . "\t" 
								. $SNAP_tagsnp_to_proxysnp{$proxy_id}{'rsquared'} . "\t" 
								. $SNAP_tagsnp_to_proxysnp{$proxy_id}{'distance'} . "\n";
				print $_, "\n\n";
			}			
		}

	}
}
close $instream;