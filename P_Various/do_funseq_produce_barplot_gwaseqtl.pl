#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;


#DRAFT 10 26/3/2015 - I am running the LD plink using GRASP (richer) and the GWAS directionality using the GWAS catalog (Grasp does not have regression parameters). Therefore I have to reintroduce both GWAS and GRASP data for the overlap. Additionally I now have the Geuvadis intersection

#not valid anymore see above
#DRAFT 6 19/1/2015 - Chris does not like the separation GWAS/GRASP and the fact that I list YRI as well. 
#Also, for the forest plot, I'm only using GWAS catalog data
#so only use GWAS catalog data, and only LD CEU data

#I want to create a bar plot for the last figure of the paper where I list one bar per category
#PICS snps
#GWAS snps
#LD GWAS
#dsQTL (LCL)
#eQTL (LCL)

#extract the snps from those annotated with funseq. Divide again by -vdr and -recur creating the following combinations
#-all
#-recur
#but also maybe:
#-vdr
#-vdr && recur

#annotation (#tag snps at max p=5E-8)
my $INPUT_GWAS = "/net/isi-scratch/giuseppe/indexes/GWAS/UCSC_gwas_catalog/gwasCatalog.b37.maxp_5E-8.minsnp5.vcf.gz";#b37
my $INPUT_GWAS_LD_CEU = "/net/isi-scratch/giuseppe/indexes/GWAS/UCSC_gwas_catalog/LD_PLINK_CEU_GWAScat.bed.gz"; #b37
my $INPUT_GRASP = "/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/GRASP2final.vcf.gz"; #remove qtls, methylation etc #b37
my $INPUT_GRASP_LD_CEU = "/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/LD_PLINK_CEU_GRASP2final.vcf.gz"; #b37

#old
#my $INPUT_GWAS = "/net/isi-scratch/giuseppe/indexes/GWAS/gwascatalog_20140218_hg19_PLUS_beechamMS_maxp_5x10-8.vcf.gz";
#my $INPUT_GRASP = "/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/GRASP2final.vcf.gz"; #remove qtls, methylation etc #b37
#my $INPUT_GWAS_LD_YRI = "/net/isi-scratch/giuseppe/indexes/GWAS/LD_PLINK_YRI_GWAScat_p_5x10-8.vcf.gz"; #b37
#my $INPUT_GRASP_LD_CEU = "/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/LD_PLINK_CEU_GRASP2final.vcf.gz"; #b37
#my $INPUT_GRASP_LD_YRI = "/net/isi-scratch/giuseppe/indexes/GWAS/GRASP/LD_PLINK_YRI_GRASP2final.vcf.gz"; #b37
my $INPUT_DSQTL = "/net/isi-scratch/giuseppe/indexes/GWAS/PRITCHARD_EQTL/DNASEqtl_degner.bed.gz";
my $INPUT_eQTL = "/net/isi-scratch/giuseppe/indexes/GWAS/PRITCHARD_EQTL/eQTL_onlylcl_v3.sorted.hg19.bed.gz";
my $INPUT_PICS = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/BROAD_PICS/BROAD_candidate_causal_snps_39immune_nonimmune_diseases_plus_enh_annot_masterfile9_hg19.vcf.gz";

#these are ALL eQTLs; the GAT test was with all eQTLs; the directionality analysis is with the subset of the best.
my $INPUT_GEUVADIS_EQTL_CEU = "/net/isi-scratch/giuseppe/indexes/GEUVADIS/EUR373.alltypes.cis.FDR5.all.rs137.GAT_hg19.bed.gz"; #hg19
my $INPUT_GEUVADIS_EQTL_YRI = "/net/isi-scratch/giuseppe/indexes/GEUVADIS/YRI89.alltypes.cis.FDR5.all.rs137.GAT_hg19.bed.gz";  #hg19
#these are the BEST
my $INPUT_GEUVADIS_EQTL_CEUb = "/net/isi-scratch/giuseppe/indexes/GEUVADIS/BEST/EUR373.alltypes.cis.FDR5.best.rs137.GAT.txt.hg19.bed.gz"; #hg19
my $INPUT_GEUVADIS_EQTL_YRIb = "/net/isi-scratch/giuseppe/indexes/GEUVADIS/BEST/YRI89.alltypes.cis.FDR5.best.rs137.GAT.txt.hg19.bed.gz";  #hg19

my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin/bedtools";


#FILTER - only consider variants marked in a VDR peak
#change this if you don't want it
#only consider variants which replicate
#test this
#also intersect with pics and GWAS and strong LD ceu and yri
my $infile;
my $FILTER_VDR_ONLY;
my $FILTER_RECUR;
GetOptions(
        'i=s'        =>\$infile,
        'vdr'   	 =>\$FILTER_VDR_ONLY,
        'recur'   	 =>\$FILTER_RECUR,
);
if(!$infile){
	print "USAGE: do_funseq_produce_barplot.pl -i=<INFILE> -vdr -recur\n";
    print "<INFILE> vcf output (no_DBRECUR) of funseq2\n";
    print "(optional)<vdr> whether to only use ASB variants in a VDR peak (default=no)\n";
    print "(optional)<recur> whether to only use ASB variants which recur (default=no)\n";
    exit 1;
}
my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
my $TEMP_variants = $directory . $basename . 'temp.vcf';
my $TEMP_variants_b37 = $directory . $basename . 'temp_b37.vcf';

#get all snps or the subset with VDR PEAK, RECUR or BOTH annotation
my %id_to_info;
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/);
	next if($_ eq '');

	
	#get info field
	my ($chr, $pos, $rsID, $info) = (split /\t/)[0,1,2,7];
	#filters-------------------------------------------------
	if($FILTER_RECUR){
		next unless($info =~ /RECUR/);
	}
	if($FILTER_VDR_ONLY){
		my $info_ncenc;
		my @info = split(";", $info);
		foreach my $item (@info){
			$info_ncenc = $item if($item =~ /^NCENC/);
		}
		next if(!$info_ncenc);
		next unless($info_ncenc =~ /TFP(.*)VDR/);
	}
	#-------------------------------------------------------
	my $id = $chr . '-' . $pos . '-' . $rsID; #recurrent is  printed once per sample in output.vcf
	$id_to_info{$id} = $info;
}

open (my $outstream,  q{>}, $TEMP_variants) or die("Unable to open $TEMP_variants : $!");
print $outstream '##fileformat=VCFv4.0' . "\n";
print $outstream '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO' . "\n";
foreach my $item (sort keys %id_to_info){
	my @item = split('-', $item);
	print $outstream $item[0] . "\t" . $item[1] . "\t" . $item[2] . "\t-\t-\t100\tPASS\tNA\n"; 
}
close $outstream;

open ($outstream,  q{>}, $TEMP_variants_b37) or die("Unable to open $TEMP_variants_b37 : $!");
print $outstream '##fileformat=VCFv4.0' . "\n";
print $outstream '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO' . "\n";
foreach my $item (sort keys %id_to_info){
	my @item = split('-', $item);
	my $chr = $item[0];
	$chr =~ s/chr(.*)/$1/;
	print $outstream $chr . "\t" . $item[1] . "\t" . $item[2] . "\t-\t-\t100\tPASS\tNA\n"; 
}
close $outstream;

my $counter = `cat $TEMP_variants | wc -l`;
$counter -= 2; #header

#intersections and report
#numbers will be number of ASB SNP with AT LEAST 1 intersection in annotation
my $out_grasp        = `$BEDTOOLS intersect  -u -a $TEMP_variants_b37 -b $INPUT_GRASP        | wc -l`;
my $out_grasp_ld_ceu = `$BEDTOOLS intersect  -u -a $TEMP_variants_b37 -b $INPUT_GRASP_LD_CEU | wc -l`;
#my $out_grasp_ld_yri = `$BEDTOOLS intersect  -u -a $TEMP_variants_b37 -b $INPUT_GRASP_LD_YRI | wc -l`;
#my $out_gwas_ld_yri  = `$BEDTOOLS intersect  -u -a $TEMP_variants_b37 -b $INPUT_GWAS_LD_YRI  | wc -l`;

my $out_gwas         = `$BEDTOOLS intersect  -u -a $TEMP_variants_b37 -b $INPUT_GWAS         | wc -l`;
my $out_gwas_ld_ceu  = `$BEDTOOLS intersect  -u -a $TEMP_variants_b37 -b $INPUT_GWAS_LD_CEU  | wc -l`;
my $out_dsqtl        = `$BEDTOOLS intersect  -u -a $TEMP_variants     -b $INPUT_DSQTL        | wc -l`;
my $out_eqtl         = `$BEDTOOLS intersect  -u -a $TEMP_variants     -b $INPUT_eQTL         | wc -l`;
my $out_pics         = `$BEDTOOLS intersect  -u -a $TEMP_variants     -b $INPUT_PICS         | wc -l`;

my $out_eqtl_geuv_ceu   = `$BEDTOOLS intersect  -u -a $TEMP_variants     -b $INPUT_GEUVADIS_EQTL_CEU         | wc -l`;
my $out_eqtl_geuv_yri   = `$BEDTOOLS intersect  -u -a $TEMP_variants     -b $INPUT_GEUVADIS_EQTL_YRI         | wc -l`;
my $out_eqtl_geuv_ceub  = `$BEDTOOLS intersect  -u -a $TEMP_variants     -b $INPUT_GEUVADIS_EQTL_CEUb         | wc -l`;
my $out_eqtl_geuv_yrib  = `$BEDTOOLS intersect  -u -a $TEMP_variants     -b $INPUT_GEUVADIS_EQTL_YRIb         | wc -l`;

#my $out_gwas         = `$BEDTOOLS intersect  -a $TEMP_variants     -b $INPUT_GWAS         | wc -l`;
#my $out_grasp        = `$BEDTOOLS intersect  -a $TEMP_variants_b37 -b $INPUT_GRASP        | wc -l`;
#my $out_gwas_ld_ceu  = `$BEDTOOLS intersect  -a $TEMP_variants_b37 -b $INPUT_GWAS_LD_CEU  | wc -l`;
#my $out_gwas_ld_yri  = `$BEDTOOLS intersect  -a $TEMP_variants_b37 -b $INPUT_GWAS_LD_YRI  | wc -l`;
#my $out_grasp_ld_ceu = `$BEDTOOLS intersect  -a $TEMP_variants_b37 -b $INPUT_GRASP_LD_CEU | wc -l`;
#my $out_grasp_ld_yri = `$BEDTOOLS intersect  -a $TEMP_variants_b37 -b $INPUT_GRASP_LD_YRI | wc -l`;
#my $out_dsqtl        = `$BEDTOOLS intersect  -a $TEMP_variants     -b $INPUT_DSQTL        | wc -l`;
#my $out_eqtl         = `$BEDTOOLS intersect  -a $TEMP_variants     -b $INPUT_eQTL         | wc -l`;
#my $out_pics         = `$BEDTOOLS intersect  -a $TEMP_variants     -b $INPUT_PICS         | wc -l`;


print "ANNOTATION\tOBSERVED\n";
print 'ALL ANNOTATED VDR VARIANTS' . "\t" . $counter . "\n";
print 'GWAS catalog' . "\t" . $out_gwas;
print 'GRASP catalog' . "\t" . $out_grasp;
print 'GWAS catalog, (strong LD, CEU)' . "\t" . $out_gwas_ld_ceu;
#print 'GWAS catalog, (strong LD, YRI)' . "\t" . $out_gwas_ld_yri;
print 'GRASP catalog (strong LD, CEU)' . "\t" . $out_grasp_ld_ceu;
#print 'GRASP catalog (strong LD, YRI)' . "\t" . $out_grasp_ld_yri;
print 'dsQTL (Degner et al. 2012)' . "\t" . $out_dsqtl;
print 'eQTL (Pritchard eQTL, LCL only)' . "\t" . $out_eqtl;
print 'eQTL (GEUVADIS, CEU)' . "\t" . $out_eqtl_geuv_ceu;
print 'eQTL (GEUVADIS, YRI)' . "\t" . $out_eqtl_geuv_yri;
print 'eQTL (GEUVADIS, CEU, best)' . "\t" . $out_eqtl_geuv_ceub;
print 'eQTL (GEUVADIS, YRI, best)' . "\t" . $out_eqtl_geuv_yrib;
print 'BROAD PICS (autoimmune causal)' . "\t" . $out_pics;

unlink $TEMP_variants;