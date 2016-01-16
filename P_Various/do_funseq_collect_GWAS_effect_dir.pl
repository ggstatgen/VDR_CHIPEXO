#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#Ram is confused by the fact that there are VDR-BVs and VDR-hcBVs
#Also, should I include cases where the OR or CIs are not reported but the rest is?

#23/3/2015 Chris wants me to discard all positions that do not have ancestral state. Maybe count them to see how many are there?

#Chris wants to see the DAF instead of the CAF from the GWAS. Adding this

#in figure s21, the counts refer to snps with p < 5E-8 and traits with at least 5 snps in total. Adjust here.

#You want to draw a summary figure including the VDR-BVs which are GWAS snps, how they affect the binding of VDR, what is the pvalue, the number of cases in the study, any close gene, maybe the LD information, and especially the DIRECTIONALITY/ODD RATIO, allele frequency, etc

#the idea is to summarize all phenotypic/genotypic/disease info for a few selected candidates and to show whether the SNP has a protective or deleterious effect for each disease.

#This post is very nice
#http://gettinggeneticsdone.blogspot.co.uk/2011/06/displaying-regression-coefficients-from.html
#and points to the following methodology:
#http://visualization.ritchielab.psu.edu/synthesis_views/plot
#paper http://visualization.ritchielab.psu.edu/synthesis_views/plot

#This script will prepare the inputs for the visualisation in the link.
#inputs: http://visualization.ritchielab.psu.edu/synthesis_views/document
#Input could be something like this:
#SNP	Chr	pos	MainGroup:pval	MainGroup:or	MainGroup:lower_ci	MainGroup:upper_ci	MainGroup:caf
#RS1801133	1	11856378	0.00000005	1.006	0.847	1.194	0.266444197			
#RS1748195	1	63049593	7E-016	0.888	0.759	1.039	0.093344876

#GWAS catalog info
#http://www.genome.gov/27529028
#For each identified SNP, we extract: chromosomal region (from UCSC Genome Browser); gene (as reported); rs number and risk allele (as reported); risk allele frequency in controls (if not available among all controls, among the control group with the largest sample size); p-value and any relevant text (e.g., subgroups where applicable); OR (or % variance explained, SD increment, or unit difference for quantitative traits), 95% CI and any relevant text (e.g., subgroups). If the p-value, OR, and 95% CI fields are not available for the combined population, we extract estimates from the population group with the largest sample size.

#you have: rsId, chr, pos, pval, OR, lower_CI, upper_ci, risk allele frequency
#to this, you will add GWAS genotype direction (protective or causative?) and phenotype direction (LOB or GOB)

#INPUT: file obtained by intersecting funseq2 output.vcf and gwascatalog.bed file obtained from UCSC GWAS bed, in turn obtained using do_GWAScatalog_UCSC_to_bed.pl

#NOTE: # is NOT A GOOD SEPARATOR, there are some in the fields it self! Change!

#INPUT DIR: /net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/GWAS_INTERSECTIONS/DIRECT_GWAS_INTERSECTION

#What to do when risk allele is neither ANC nor ALT? 
#Reverse complement, and assume value given on wrong strand

#my $INPUT_VARIANTS = '/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/ALLELESEQ_20samples_merged.vcf.gz'; #b37
my $INPUT_VARIANTS = "/net/isi-mirror/1000_Genomes/ALL.phase1_release_v3.20101123.snps_nogts.vcf.gz";
my $INPUT_ALLELESEQ_RAW = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/_RAW/interestingHets_vdrchipexo.txt"; #b37
my $INPUT_GWAScat = "/net/isi-scratch/giuseppe/indexes/GWAS_Gwascat/UCSC_gwas_catalog/gwasCatalog.b37.BeechamMS.bed"; #b37

my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";

#loss of binding, gain of binding
my $LOB = 'lob'; #ancestral -> alternative reduces binding
my $GOB = 'gob'; #ancestral -> alternative increase binding
#loss of risk, gain of risk
my $LOR = 'lor'; #ancestral -> alternative reduces risk  (i.e. risk allele is ancestral a.)
my $GOR = 'gor'; #ancestral -> alternative increases risk  (i.e. risk allele is alternative a.)

#LOB + LOR -> transition from ancestral to alternative cause:  1) decrease of VDR binding 2) decrease of disease risk - VDR-BV is DELETERIOUS
#GOB & GOR -> transition from ancestral to alternative cause:  1) increase of VDR binding 2) increase of disease risk - VDR-BV is DELETERIOUS

#LOB & GOR -> transition from ancestral to alternative cause:  1) decrease of VDR binding 2) increase of disease risk - VDR is PROTECTIVE
#GOB & LOR -> transition from ancestral to alternative cause:  1) increase of VDR binding 2) decrease of disease risk - VDR is PROTECTIVE


my $infile;
my $output;
my $subset_bed; #if you only want to collect data for positions in an external bed file (eg CLASS I sites)
my $title;
GetOptions(
        'i=s'      =>\$infile,
        'subset=s' =>\$subset_bed,
        'title=s'  =>\$title
);
#$infile = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/GWAS_INTERSECTIONS/DIRECT_GWAS_INTERSECTION/Output_hg19_allsamples_INTERSECT_gwasCatalog.hg19.tsv"; #hg19
#$subset_bed ="/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_motif_scan_pscanchip/Pscanchip_occurrences_VDRRXR_jaspar_classI_peakintervals.bed"; #hg19

#allow to select items in class I sites only
if(!$infile){
     print "USAGE: do_funseq_collect_GWAS_effect_dir.pl -i=<INFILE> -subset=<SUBSET_BED> -title=<TITLE>\n";
     print "<INFILE> [b37] vcf.gz with the variants you want to test (check FIG6_*.vcf.gz files)\n";
     print "(optional)<SUBSET_BED> only collect data for VDR-BVs intersecting with this bed file(default=no)\n";
     print "(optional)<TITLE> when SUBSET_BED is set, this may contain a suffix for the output file\n";
     print "Gwas variants will be from $INPUT_GWAScat. Change if needed\n";
     exit 1;
}
my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
if($subset_bed){
	$output  = $directory . $basename . '_gwas_effect'  . $title .  '.Rdata';
}else{
	$output  = $directory . $basename . '_gwas_effect.Rdata';
}
my $input_intersect_gwas = $directory . 'TEMP_gwas_effect_' . $basename . 'data';


#################
#1 get intersection
##################
system "$BEDTOOLS intersect -wb -a $infile -b $INPUT_GWAScat > $input_intersect_gwas";

################
#2 if there is a bed file with peaks (eg class I peaks) save their coordinates in hash
################
my %subset_bed;
if($subset_bed){
	open (my $instream,  q{<}, $subset_bed) or die("Unable to open $subset_bed : $!");
	while(<$instream>){
		chomp;
		my ($chr, $start, $stop) = (split /\t/)[0,1,2];
		my $interval = $start . '-' . $stop;
		$subset_bed{$chr}{$interval} = 1;
	}
	close $instream;
}
#####################
#3 -build hash to map chr-pos to ref, ancestral and alternate
#####################
#chr6    32590953        rs111410428;rs9271588
#does not have an ancestral, but I found an ancestral on ensembl variation
#T/C ancestral T

my %variants_1kg;
tie *FILE,   'IO::Zlib', $INPUT_VARIANTS, "rb";
while (<FILE>)	{ 
	chomp;
	next if($_ eq '');
	next if($_ =~ /^#/);
	#my $ALTF_AFR; 
	my $AF_EUR; my $AF_GLOB; #ALTERNATE allele frequencies
	my $FLAG; # set to one if the ancestral is the alternate
	my $ANC;	
		
	my @fields = split("\t", $_);
	#skip indels
	next if(length($fields[3]) > 1);
	next if(length($fields[4]) > 1);
	
	my $key = $fields[0] . '-' . $fields[1];
	my $ref = uc($fields[3]); 
	my $alt = uc($fields[4]);
	my $info = $fields[7];	
	
	#get ancestral
	my @info = split(";", $fields[7]);	
	
#	##chr6    32590953        rs111410428;rs9271588
#	if($key eq '6-32590953'){
#		$info[0] = 'AA=T';
#	}
	
	#get ancestral allele info and allele frequencies=====================
	foreach my $item (@info){
		if($item =~ /^AA=(.*)/){
			$ANC = $1;
		}
		if($item =~ /^EUR_AF=(.+)/){
			$AF_EUR = $1;	
		}
		if($item =~ /^AF=(.+)/){
			$AF_GLOB = $1;
		}	
	}	
	
	#if no AA info, or ambiguities, move on
	next unless($ANC);
	if(uc($ANC) eq uc($alt)){
		$FLAG = 1;	
		$variants_1kg{$key}{ANC} = $alt;
		$variants_1kg{$key}{DER} = $ref;		
	}elsif(uc($ANC) eq uc($ref)){
		$variants_1kg{$key}{ANC} = $ref;
		$variants_1kg{$key}{DER} = $alt;					
	}elsif( ($ANC eq '') or ($ANC eq 'N') or ($ANC eq '.') or ($ANC eq '-') ){
		next;
	}else{ 
		next;			
	}	

	if($AF_EUR){
		if($FLAG){
			if($AF_EUR == 1){
				$variants_1kg{$key}{DAFEUR} = '0.0';				
			}else{
				$variants_1kg{$key}{DAFEUR} = (1 - $AF_EUR);				
			}
		}else{
			$variants_1kg{$key}{DAFEUR} = $AF_EUR;
		}
	}
	if($AF_GLOB){
		if($FLAG){
			if($AF_GLOB == 1){
				$variants_1kg{$key}{DAFGLOB} = '0.0';				
			}else{
				$variants_1kg{$key}{DAFGLOB} = (1 - $AF_GLOB);				
			}
		}else{
			$variants_1kg{$key}{DAFGLOB} = $AF_GLOB;
		}
	}	
}
close FILE;

###################
#4 collect VDR-BV phenotypes, get also pvalues
###################
#warning about non matching alleles
#6-32345443: (ancestral ref/alt and mat/pat don't match: (G, A) and (G, T)
#chr-position => samplename => alleles / read counts info
my %position2sample2readdepth;
open (my $instream,  q{<}, $INPUT_ALLELESEQ_RAW) or die("Unable to open $INPUT_ALLELESEQ_RAW : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^sample/); #header
	next if($_ eq '');
	
	#format
	#sample	chrm	snppos	ref	mat_gtyp	pat_gtyp	c_gtyp	phase	mat_all	pat_all	cA	cC	cG	cT	winning	SymCls	SymPval	BindingSite	cnv
	#NA06986	1	1080920	G	S	W	R	PHASED	G	A	2	0	7	0	M	Sym	0.1796875	1	1.0
	my ($sample_id, $chr,$snppos, $ref, $m_allele, $p_allele, $cA, $cC, $cG, $cT, $winner, $symcls, $asb_pval) = (split /\t/)[0,1,2,3,8,9,10,11,12,13,14,15,16];
	
	next if(!$chr);
	next if(!$snppos);
	next unless ($symcls =~ /Asym/);
	my $coordinate_id = $chr . '-' . $snppos;
	
	my $allele_ancestral = $variants_1kg{$coordinate_id}{ANC};
	my $allele_derived = $variants_1kg{$coordinate_id}{DER}; 
	next unless($allele_ancestral);
	next unless($allele_derived);

	if( ($m_allele eq 'None') or ($p_allele eq 'None')  ){
	}else{
		if( ($allele_ancestral eq $m_allele)  && ($allele_derived eq $p_allele) ){		
		}elsif( ($allele_ancestral eq $p_allele)  && ($allele_derived eq $m_allele)  ){	
		}else{
			print STDERR "$coordinate_id: (ancestral ref/alt and mat/pat don't match: ($allele_ancestral, $allele_derived) and ($m_allele, $p_allele)\n";
			print STDERR "Most probably the 1kg info is from an inputed data point. Saving paternal or maternal info\n";
			next;
		}		
	}
	my $read_data = join(",", $allele_ancestral,$allele_derived, $cA, $cC, $cG, $cT, $asb_pval);
	$position2sample2readdepth{$coordinate_id}{$sample_id} = $read_data;
}
close $instream;

#INPUT FORMAT
#[0,7] funseq vcf
#[8-11] gwas bed
#0	chr
#1	VDR-BV position
#2	VDR-BV rsID
#3	VDR-BV ref
#4	VDR-BV alt
#5 	100
#6	PASS
#7	INFO (from here get motifbr and link to reads to get directionality)
#8	chr
#9	start
#10	stop
#11	INFO GWAS. This is # separated. Split and process to obtain pval, OR, lower_CI, upper_CI, allele freq

#FIELD 11, subfields
#11.0 rsID #test it's same as 2
#11.1 PUBMED ID
#11.2 AUTHOR
#11.3 date
#11.4 JOURNAL
#11.5 title
#11.6 TRAIT
#11.7 initial sample size
#11.8 rep sample size
#11.9 region
#11.10 gene
#11.11	strongest snp risk allele
#11.12	risk allele frequency - can be NR, not reported
#11.13	p-value
#11.14	?
#11.15	OR
#11.16	95 CI text
#11.17	Platform
#11.18	CNV

#SNP	Chr	pos	MainGroup:pval	MainGroup:or	MainGroup:lower_ci	MainGroup:upper_ci	MainGroup:caf

open (my $outstream,  q{>}, $output) or die("Unable to open $output : $!");
print $outstream "SNP\tChr\tpos\tMainGroup:pval\tMainGroup:or\tMainGroup:lower_ci\tMainGroup:upper_ci\tMainGroup:caf\tMainGroup:CEU_DAF\tCI_DATA\tBINDING_MAGNITUDE\tLCL_SAMPLE\tASB_PVAL\tBINDING_DIR\tRISK_DIR\tINIT_SAMPLE\tREP_SAMPLE\tEFFECT\tTRAIT_DISEASE\tAUTHOR\tSENTINEL_GWAS_SNP\tRC\n"; 
open ($instream,  q{<}, $input_intersect_gwas) or die("Unable to open $input_intersect_gwas : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/);
	next if($_ eq '');
	my ($chr_vdrbv,$pos_vdrbv,$rsid_vdrbv,$ref_vdrbv,$alt_vdrbv,$unused1,$unused2,$info_vdrbv,$chr_gwas,$start_gwas,$stop_gwas,$info_gwas) = split("\t",$_);
	if($subset_bed){ next unless(defined check_coords_in_bed($chr_vdrbv, $pos_vdrbv)); }
	
	my $chr_b37 = $chr_vdrbv; 
	$chr_b37 =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19
	my $genomic_coord = $chr_b37 . '-' . $pos_vdrbv;
	if(!$position2sample2readdepth{$genomic_coord}){
		print STDERR "$rsid_vdrbv - $genomic_coord: no ancestral/derived info available for this VDR-BV. Skipping..\n";
		next;
	}		
	my $phenotype_direction; my $gwas_risk_direction; my $EFFECT; my $gwas_risk_is_reverse_strand; my $sample_id;
	#==================================
	#get direction and magnitude of phenotypic effect
	#==================================
	#you need to query %position2sample2readdepth so you'll need coordinates of the snp and sample id
	my ($phenotype_fold_change, $ASB_PVAL_THISSAMPLE) = get_phenotype_fc_direction($chr_vdrbv,$pos_vdrbv,$rsid_vdrbv,$info_vdrbv,\$sample_id);	
	if(!$phenotype_fold_change){
		print STDERR "$rsid_vdrbv - Warning: fold change value is empty or undefined. Skipping..\n";
		next;
	}
	if($phenotype_fold_change > 1){ #(rc_anc + 1) > (rc_der + 1) 
		$phenotype_direction = $LOB;
	}elsif($phenotype_fold_change < 1){ #(rc_anc + 1) < (rc_der + 1)
		$phenotype_direction = $GOB;
	}else{ #(rc_anc + 1) = (rc_alt + 1)
		print STDERR "$rsid_vdrbv - Warning: fold change is $phenotype_fold_change for VDR-BV. Skipping..\n";
		next;
	}
	my $BINDING_MAGNITUDE = abs($phenotype_fold_change);
	my $DAF_CEU = '-';
	if($variants_1kg{$genomic_coord}{DAFEUR}){
		$DAF_CEU = $variants_1kg{$genomic_coord}{DAFEUR};
	}else{
		print STDERR "$rsid_vdrbv - Warning: CEU DAF is not available for VDR-BV\n";
	}
	#==================================
	#risk allele: is it the ancestral or the alternative?
	#==================================
	my @gwas_info = split('#',$info_gwas);
	#sometimes in $rsid_vdrbv there is more than one rsid, test both
	if($rsid_vdrbv =~ /;/){
		my $found; my @rs_ids = split(";",$rsid_vdrbv); 
		foreach my $item (@rs_ids){
			$found = 1 if($item eq $gwas_info[0]);
		}
		unless($found){
			print STDERR "$rsid_vdrbv - Warning: the VDR-BV rsID(s) and the GWAS rsID: $gwas_info[0] don't match. Skipping..\n";
			next;
		}		
	}else{
		unless($rsid_vdrbv eq $gwas_info[0]){
			print STDERR "$rsid_vdrbv - Warning: the VDR-BV rsID(s) and the GWAS rsID: $gwas_info[0] don't match. Skipping..\n";
			next;
		}		
	}
	my $gwas_author = $gwas_info[2];
	my $gwas_trait = $gwas_info[6];
	my $gwas_init_sample = $gwas_info[7];
	my $gwas_rep_sample = $gwas_info[8];
	my $gwas_risk_allele = $gwas_info[11]; #postprocess - typical: rs8080944-A or rs8080944-?
	my $gwas_risk_allele_freq = $gwas_info[12]; #postprocess - typical 0.xxx , NR etc
	my $gwas_pval = $gwas_info[13];
	my $gwas_OR = $gwas_info[15]; #postprocess - typical: 1.45, NR etc
	my $gwas_CI = $gwas_info[16]; #postprocess - typical: [0.xx-1.33](text), NR
	
	print STDERR "$rsid_vdrbv - Warning: OR is not reported\n" if($gwas_OR eq 'NR');
	$gwas_risk_direction = get_gwas_risk_direction($chr_vdrbv,$pos_vdrbv, $rsid_vdrbv, $gwas_risk_allele,\$gwas_risk_is_reverse_strand);
	if(!$gwas_risk_direction){
		print STDERR "$rsid_vdrbv - Warning: unable to determine risk allele direction. Skipping..\n";
		next;
	}
	if(  ( ($phenotype_direction eq $LOB) && ($gwas_risk_direction eq $LOR) ) || ( ($phenotype_direction eq $GOB) && ($gwas_risk_direction eq $GOR) ) ){
		$EFFECT = 'deleterious';
	}elsif(  ( ($phenotype_direction eq $LOB) && ($gwas_risk_direction eq $GOR) ) || ( ($phenotype_direction eq $GOB) && ($gwas_risk_direction eq $LOR) ) ){
		$EFFECT = 'protective';		
	}else{
		print STDERR "$rsid_vdrbv: directionality problem. It shouldn't be here. Skipping..\n";
		next;
	} 	
	#split CIs
	my $low_CI = ''; my $hi_CI = '';
	get_CI_intervals($gwas_CI, \$low_CI, \$hi_CI);

#print $outstream "SNP\tChr\tpos\tMainGroup:pval\tMainGroup:or\tMainGroup:lower_ci\tMainGroup:upper_ci\tMainGroup:caf\tMainGroup:CEU_DAF\tCI_DATA\tBINDING_MAGNITUDE\tBINDING_DIR\tRISK_DIR\tINIT_SAMPLE\tREP_SAMPLE\tEFFECT\tTRAIT_DISEASE\tAUTHOR\tRC\n";
	
	print $outstream  	$rsid_vdrbv   . "\t" .
						$chr_vdrbv . "\t" .
						$pos_vdrbv . "\t" . 
						$gwas_pval . "\t" . 
						$gwas_OR   . "\t" .
						$low_CI  . "\t" .
						$hi_CI   .   "\t" .
						$gwas_risk_allele_freq . "\t" .
						$DAF_CEU   . "\t" .
						$gwas_CI   . "\t" .
						$BINDING_MAGNITUDE . "\t" .
						$sample_id         .  "\t" . 
						$ASB_PVAL_THISSAMPLE . "\t" . 												
						$phenotype_direction . "\t" .
						$gwas_risk_direction . "\t" .
						$gwas_init_sample . "\t" . 
						$gwas_rep_sample . "\t" . 
						$EFFECT . "\t" . 
						$gwas_trait . "\t" . 
						$gwas_author . "\t" . 
						$gwas_risk_is_reverse_strand . "\n";
}
close $instream;
close $outstream;
unlink $input_intersect_gwas;

#############
#SUBROUTINES
#############
#------------------------------------------------------------------------
sub get_gwas_risk_direction{
	my ($chr, $pos, $ID, $risk_allele_exp, $is_reverse_strand) = @_;
	$$is_reverse_strand = '';
	
	#clean up risk allele
	my $risk_allele; 
	if($risk_allele_exp =~ /rs(\d+)-([AGTC\?])/){
		$risk_allele = $2;
	}elsif($risk_allele_exp =~ /rs(\d+)-([agct\?])/){
		$risk_allele = uc($2);
	}else{
		print STDERR "$ID - get_gwas_direction(): the R.A. is in an expression I don't understand: $risk_allele_exp. Skipping..\n";
		return undef;
	}
	return undef if($risk_allele eq '?');

	$chr =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19
	my $coord_id = $chr . '-' . $pos;
	my $anc_ref = $variants_1kg{$coord_id}{ANC};
	my $anc_alt = $variants_1kg{$coord_id}{DER};
	#next unless($anc_ref);
	#next unless($anc_alt);	
	
	
	#is the risk allele the ancestral or the derived?
	if($risk_allele eq $anc_ref){ #risk allele is the ancestral reference. ancestral -> derived = lor
		#print STDERR "$ID - get_gwas_direction(): the R.A. is the ancestral: $risk_allele\n";
		return $LOR;
	}elsif($risk_allele eq $anc_alt){#risk allele is the alternative. ancestral -> derived = gor
		#print STDERR "$ID - get_gwas_direction(): the R.A. is the derived: $risk_allele\n";	
		return $GOR;
	}else{
		print STDERR "$ID - get_gwas_direction(): the R.A is neither the ancestral nor the derived: (R.A., anc, der) = ($risk_allele, $anc_ref, $anc_alt)\n";
		print STDERR "$ID - get_gwas_direction(): testing the reverse complement to see if the RA was reported on the other strand..\n";	
		if($risk_allele eq 'A'){
			$risk_allele = 'T';
		}elsif($risk_allele eq 'T'){
			$risk_allele = 'A';
		}elsif($risk_allele eq 'G'){
			$risk_allele = 'C';
		}elsif($risk_allele eq 'C'){
			$risk_allele = 'G';
		}else{
			print STDERR "$ID - get_gwas_direction(): risk allele not recognised: $risk_allele\n. Aborting..\n"; # do this check above.
			return undef;
		}
		if($risk_allele eq $anc_ref){
			print STDERR "$ID - get_gwas_direction(): OK the R.A. is now the ancestral: $risk_allele\n";
			$$is_reverse_strand = 'Y';
			return $LOR;			
		}elsif($risk_allele eq $anc_alt){
			print STDERR "$ID - get_gwas_direction(): OK the R.A. is now the derived: $risk_allele\n";
			$$is_reverse_strand = 'Y';
			return $GOR;
		}else{
			print STDERR "$ID - get_gwas_direction(): the R.A is neither the ancestral nor the derived even after reversing: (R.A., anc, der) = ($risk_allele, $anc_ref, $anc_alt). Giving up..\n";
			return undef;
		}		
	}	
}

#------------------------------------------------------------------------
sub get_phenotype_fc_direction{
	my ($chr, $pos, $rs_id, $info, $sample_id) = @_;
	
	#$chr =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19
	my $genomic_coord = $chr . '-' . $pos;
	my @info = split(";", $info);
	#get sample ID
	#typical:SAMPLE=interestingHets_NA07029_EBLfiltered_hg19
	if($info[0] =~ /SAMPLE=interestingHets\_(NA\d{5})\_EBLfiltered/){
		$$sample_id = $1;
	}else{
		print STDERR "Error: I cannot find the requested sample: $info[0] in the input information. Aborting\..n";
		exit -1;
	}
	
	if(!$position2sample2readdepth{$genomic_coord}{$$sample_id}){
		print STDERR "$rs_id - Warning no read depth data found for VDR-BV at position: $genomic_coord, sample: $$sample_id\n";
		return undef;
	}
	my ($fc, $pval) = process_read_data($position2sample2readdepth{$genomic_coord}{$$sample_id}, $$sample_id, $rs_id);
	if(!$fc || $fc eq ''){
		print STDERR "Warning fold change value $fc is empty or undefined for VDR-BV: $rs_id. Skipping..\n";
		return undef;
	}else{
		return ($fc,$pval);
	}
}


#------------------------------------------------------------------------
#this splits a read_data string of the kind "c,a,,0,6,0,0"
#into a phenotype disruption score FC = [ (alt_count + 1) / (anc_count + 1) ] (turned into logs in R)
#I do alt/anc so that LOB points are negative on the graph, while GOB are positive
#if the motif hitting snp with the alternative allele has lowered the coverage at the position, this subroutine should output a negative value

#input meaning:
#read data string format: ancestral_ref, ancestral_alt, $cA,$cC,$cG,$cT;
sub process_read_data{
	my ($read_data_string, $lcl_sample, $id) = @_;
	
	my ($allele_anc,$allele_der,$cA,$cC,$cG,$cT, $thispval) = split(",", $read_data_string);
	my $anc_count; my $der_count;
	my $anc_base = $allele_anc;
	my $der_base = $allele_der;

	if( ($cA eq 0) && ($cC eq 0) && ($cG eq 0) && ($cT eq 0) ){
		print STDERR "process_read_data(): error asb snp: $id. (LCL sample: $lcl_sample) - 0 read coverage for all possible nucleotides at this event ($cA,$cC,$cG,$cT). Aborting..\n";
		exit -1;
	}

	#get reference read coverage
	if($anc_base eq 'A'){
		$anc_count = $cA;
	}elsif($anc_base eq 'C'){
		$anc_count = $cC;
	}elsif($anc_base eq 'G'){
		$anc_count = $cG;
	}elsif($anc_base eq 'T'){
		$anc_count = $cT;
	}else{
		print STDERR "process_read_data(): error asb snp: $id. (LCL sample: $lcl_sample) - reference base not recognised: $anc_base. Aborting..\n";
		exit -1;
	}
	
	if($der_base eq 'A'){
		$der_count = $cA;
	}elsif($der_base eq 'C'){
		$der_count = $cC;
	}elsif($der_base eq 'G'){
		$der_count = $cG;
	}elsif($der_base eq 'T'){
		$der_count = $cT;
	}else{
		print STDERR "process_read_data(): error asb snp: $id. (LCL sample: $lcl_sample) - alternate base not recognised: $der_base. Aborting..\n";
		exit -1;
	}
	my $fc =  ( $anc_count + 1 ) / ( $der_count + 1 );
	return ($fc, $thispval);
}

#------------------------------------------------------------------------
sub check_coords_in_bed{
	my ($this_br_chr, $this_br_pos) = @_;
	
	foreach my $interval (keys %{ $subset_bed{$this_br_chr} }){
		my ($peak_start, $peak_end) = split('-', $interval);
		if( ($this_br_pos <= $peak_end)  &&  ($this_br_pos >= $peak_start) ){
			return 1;
		}
	}
	return undef;
}
#----------------------------------------------------------------------------

#[1.23-1.45]
#[1.71-4.31] unit increase
#NR
sub get_CI_intervals{
	my ($ci_data, $low_ci, $high_ci) = @_;
	return if($ci_data =~ /NR/);
	if($ci_data =~ /\[(.+)-(.+)\](.*)/){
		$$low_ci = $1;
		$$high_ci = $2;
	}else{
		return;
	}
}