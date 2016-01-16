#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;
use List::Util qw(min max sum reduce);

##update 26-8-2015
#if REP, you need to remove all instances where direction of binding differs across samples

###update
#you need to be able to select variants based on ALLELESEQ threshold

###FIGURE 6A/6B
#Script to get independent sets of VDR-BVs where I select only one for each D' LD block based on the following:
#the VDR-BV causes the largest magnitude of binding affinity change (median across samples of difference of binding affinity change)
#if more than one VDR-BV has same magnitude, collect the one with largest amount of maximum mapped reads.

#The script also restricts the background to independent variants [options: 1) variant CHOSEN ABOVE OR 2) random variant in block]

#before carrying out the LD clumping, the script removes any foreground and background variables in the MHC defined by:
#/net/isi-scratch/giuseppe/VDR/MHC_b37.bed
#6	29540169	33215544	MHC_REGION

#THOUGHTS:
#There are many ways you could prioritise further these VDR-BVs until there is only one of them in the LD block. Please modify the script if you find better ways.
#eg using funseq annotation. Only get the one that is in an ENCODE peak? ENCODE motif? Specific DAF value?

#10/6/2016
#Both Ram and the audience at GDS talk in Bristol asked if I had considered LD between the SNPs in the foreground.

#FOREGROUND: it will contain the BEST VDR-hcBV according to the def above.
#BACKGROUND: for each LD block, leave only the background variant corresponding to the one chosen above, or alternatively a random one.

my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";

my $in_foreground; 
my $in_background;

#the following two are needed if a p-value on the ASB is required. Build a hash out of the following two. 
my $INPUT_ALLELESEQ_RAW = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/_RAW/interestingHets_vdrchipexo.txt"; #b37
#The raw QTL calls are needed because they cannot be threshold like the ASB calls
my $INPUT_VDRQTL = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/VDR_QTL_FINAL/VDRQTL_SNPTEST_BIMBAM_o3_BFthrs_1.0_hg19_3.vcf"; #hg19, convert

my $INPUT_VARIANTS = "/net/isi-mirror/1000_Genomes/ALL.phase1_release_v3.20101123.snps_nogts.vcf.gz"; #b37 used to get a vcf out of the final bg bed
my $FOREGROUND_SET;
my $BACKGROUND_SET;
my $INPUT_LD;
my $ld_block_type;
my $PVAL_THRS;
my $pval_thrs_label = '_pthrs_';

#mhc exclusion-------------
my $MHC_chr = '6';
my $MHC_start = '29540169';
my $MHC_end = '33215544';
#mhc exclusion-------------

GetOptions(
        'fg=s'      =>\$in_foreground,
        'bg=s'      =>\$in_background,
        'fgs=s'     =>\$FOREGROUND_SET,
        'bgs=s'     =>\$BACKGROUND_SET,
        't=s'       =>\$PVAL_THRS,
        'ld=s'      =>\$ld_block_type
);
#$in_foreground = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/Output_noDBRECUR.vcf"; #hg19
#$in_background = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/d_GAT_BACKGROUNDS/GAT_BACKGROUND_SIMASYM_FINAL2.vcf.gz";
#$FOREGROUND_SET = "REP";
#$BACKGROUND_SET = "NR";
#$ld_block_type = 'D';

my $USAGE = "USAGE: do_funseq_VET_01_LDclump.pl -fg=<FG> -bg=<BG> -t=<PVAL_THRS> -fgs=<ALL|HCBV|REP|ENH> -bgs=<R|NR> -ld=<D|R>\n" .
			"<FG> [hg19] vcf output of funseq2 with noRECUR\n" .
			"<BG> [b37] vcf.gz file containing your choice of variants for GAT background\n" .
			"<PVAL_THRS> (optional) [0-0.9] if set, filter out of FG variants with ASB p-val >PVAL_THRS\n" . 
			"fgs: foreground set - one of ALL (all VDR-BVs), HCBV (RECUR), REP (only replicating), ENH (in enhancer)\n" . 
			"bgs: for each LD block containing a FG variant, pick the FG variant [NR] or a random one [R] ?\n" .  
			"ld: whether to use plink D'[D] or Mig Rsquare >0.8 [R] definition of LD block\n";

#CONVERT ALL TO b37
if(!$in_foreground){
	print $USAGE;
    exit 1;
}
if(!$in_background){
	print $USAGE;
    exit 1;
}
if(!$FOREGROUND_SET){
	print $USAGE;
    exit 1;	
}
if(!$BACKGROUND_SET){
	print $USAGE;
    exit 1;		
}
unless($FOREGROUND_SET eq 'ALL' || $FOREGROUND_SET eq 'HCBV' || $FOREGROUND_SET eq 'REP' || $FOREGROUND_SET eq 'ENH'){
	print $USAGE;
    exit 1;
}
unless($BACKGROUND_SET eq 'R' || $BACKGROUND_SET eq 'NR'){
	print $USAGE;
    exit 1;	
}
if(!$ld_block_type){
	print $USAGE;
    exit 1;		
}
unless($ld_block_type eq 'D' || $ld_block_type eq 'R'){
	print $USAGE;
    exit 1;		
}

my $ld_label;
if($ld_block_type eq 'D'){
	$INPUT_LD = "/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_CEU_LD_b37_noindels_MAF0.01.bed.gz"; #general LD blocks #b37
	$ld_label = '_LD_Dprime';	
}else{
	$INPUT_LD = "/net/isi-mirror/1000_Genomes/POPULATION_PANELS/LDEXPL_HAPLOTYPE_BLOCKS_r2/LD_r0.8_MIG_EUR_noindels_maf0.01_hwe0.001_b37.bed.gz";
	$ld_label = '_LD_r2_0.8';	
}

my($basename, $directory) = fileparse($in_foreground);
$basename =~ s/(.*)\..*/$1/;
########
#outputs
########
my $out_temp_fg_list    = $directory .'LDCLUMP_init_' .  $FOREGROUND_SET . '_' . $basename  . '_fg_list.bed'; 
my $out_temp_bg_list    = $directory .'LDCLUMP_init_' . $basename  . '_bg_list.bed';
#The above will contain unique FG positions and unique BG variable position. Used by bedtools.

my $out_temp_fg_buffer =  $directory . $basename . '_temp_fg_buffer.bed';
my $out_temp_bg_buffer =  $directory . $basename . '_temp_bg_buffer.bed'; 

my $out_temp_fg2_buffer =  $directory . $basename . '_temp_fg2_buffer.bed';
my $out_temp_bg2_buffer =  $directory . $basename . '_temp_bg2_buffer.bed'; 

my $output_fg;
my $output_fg_vcf; #needed for FIGURE 6, where we intersect vdr-bvs with GWAScatalog SNPs
my $output_fg_vcf_clean; #needed to produce a vcf of the final clean VDR-rBV set mentioned in the paper


#it will need to contain all variables BEFORE ld clumping and after filtering by pval or reproducibility
if($PVAL_THRS){
	$output_fg           =  $directory . 'LDCLUMP_' . $FOREGROUND_SET . '_' . 'FG_VDRBVs_' .  $basename  . $ld_label . $pval_thrs_label . $PVAL_THRS .  '.bed';	
	$output_fg_vcf       =  $directory . 'FIG6_'       .$basename .  '_' . $FOREGROUND_SET  . $pval_thrs_label . $PVAL_THRS .  '.vcf';
	$output_fg_vcf_clean =  $directory . 'SUPPL_DATA_' .$basename .  '_' . $FOREGROUND_SET  . $pval_thrs_label . $PVAL_THRS .  '.vcf';
}else{
	$output_fg           =  $directory . 'LDCLUMP_' . $FOREGROUND_SET . '_' . 'FG_VDRBVs_' .  $basename  . $ld_label .  '.bed';
	$output_fg_vcf       =  $directory . 'FIG6_'       . $basename .  '_' . $FOREGROUND_SET   .  '.vcf';
	$output_fg_vcf_clean =  $directory . 'SUPPL_DATA_' . $basename .  '_' . $FOREGROUND_SET   .  '.vcf';
}
my $output_bg                =  $directory . 'LDCLUMP_BG_VDRBVs_' . $basename . '_' . $BACKGROUND_SET        . $ld_label .  '.bed';
my $output_bg_vcf            =  $directory . 'LDCLUMP_BG_VDRBVs_' . $basename . '_' . $BACKGROUND_SET        . $ld_label .  '.vcf.gz';
#the output should be vcf as well, because you need the allele frequency of the variables

###########
# GLOBAL hashes
###########
my %fg_vcf;
my %fg_vcf_rep;
my %data_fg;
my %data_bg;
my %data_fg_nold;
my %data_bg_nold;
my %data_ld2fgvar; 
my %data_ld2bgvar;
my %position2sample2readdepth;
my %candidate_foreground_vars;
my %candidate_background_vars;

my %accepted_fg_vars; #if #$PVAL_THR is set, the good fg variants will be here

#############
#1.0 If ASB p-value thresholding is required, build a map of acceptable VDR-ASB PLUS VDR-QTLs to follow up in the next steps
#############
if($PVAL_THRS){
	print "*THRESHOLDING ALLELESEQ PVALS: ON*\n";
	print "Thresholding foreground variants on Alleleseq p-val: $PVAL_THRS\n";
	my $counts = 0;
	
	#get the asym
	open (my $instream,  q{<}, $INPUT_ALLELESEQ_RAW) or die("Unable to open $INPUT_ALLELESEQ_RAW : $!");
	foreach(<$instream>){
		chomp;
		next if($_ =~ /^sample/); #header
		next if($_ eq '');
		
		my($chr, $pos, $symcls, $sympval) = (split /\t/)[1,2,15,16];
		next if(!$chr);
		next if(!$pos);
		next unless ($symcls =~ /Asym/);
		#thresholding
		next unless($sympval <= $PVAL_THRS);
		
		my $coord =  $chr . '-' . $pos; #b37
		$accepted_fg_vars{$coord}++;	
	}
	close $instream;
	$counts = scalar keys %accepted_fg_vars;
	print "After thresholding, $counts unique VDR-ASB positions retained from alleleseq set\n";
	
	#get all the VDR-QTL
	my %qtl_only;
	open ($instream,  q{<}, $INPUT_VDRQTL) or die("Unable to open $INPUT_VDRQTL : $!");
	foreach(<$instream>){
		chomp;
		next if($_ =~ /^\#/); #vcf header
		next if($_ eq '');
		
		my($chr, $pos) = (split /\t/)[0,1];
		$chr =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19
		my $coord =  $chr . '-' . $pos;
		$accepted_fg_vars{$coord}++;
		$qtl_only{$coord}++;			
	}	
	close $instream;	
	$counts = scalar keys %qtl_only;
	print "$counts unique VDR-QTL positions retained\n";
	$counts = scalar keys %accepted_fg_vars;
	print "$counts unique VDR-BV positions retained.\n";
}

################
#1.1 get a bed of the variants in the foreground set from the output of funseq
#get rid of the variants in the MHC region
#report how many you got rid of
#report how many are left 
################
open (my $instream,   q{<}, $in_foreground) or die("Unable to open $in_foreground : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');	
	next if($_ =~ /^\#/);
	my $info_ncenc; my $info_motifbr; my $info_recur;
	
	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = (split /\t/)[0,1,2,3,4,5,6,7]; #chr is hg19
	$chr =~ s/chr(.*)/$1/; #this needs to be b37		
	my $fg_var = $chr . '-' . $pos;	
	my $newline = join("\t",$chr,$pos,$id,$ref,$alt,$qual,$filter,$info);
	#my $newline = join("\t",$chr,$pos,$id,$ref,$alt,$qual,$filter,'SEE_FULL_FUNSEQ_OUTPUT');
	
	if($PVAL_THRS){
		next unless($accepted_fg_vars{$fg_var});
	}
	
	my @info = split(";", $info);
	foreach my $item (@info){
		$info_ncenc = $item if($item =~ /^NCENC/); 
		$info_motifbr = $item if($item =~ /^MOTIFBR/); #any motif
		$info_recur = $item  if($item =~ /^RECUR/);			
	}
	#FLAGS
	if($FOREGROUND_SET eq 'ALL'){
		$data_fg{$fg_var} = 1;
		$fg_vcf{$newline} = 1;
		next;
	}elsif($FOREGROUND_SET eq 'REP'){
		$data_fg{$fg_var}++;
		$fg_vcf{$newline} = 1;
		next;	
	}elsif($FOREGROUND_SET eq 'HCBV'){
		if($info_recur){
			$data_fg{$fg_var} = 1;	
			$fg_vcf{$newline} = 1;		
		}else{
			next;
		}
	}elsif($FOREGROUND_SET eq 'ENH'){
		next unless($info_ncenc);
		if($info_ncenc =~ /Enhancer/){
			$data_fg{$fg_var} = 1;
			$fg_vcf{$newline} = 1;				
		}else{
			next;				
		}
	}else{
		print STDERR "Error: unrecognised: $FOREGROUND_SET. Aborting..\n";
		exit -1;
	}
}
close $instream;


my %discordant_vdrrbv_pos;
my $vcf_pos_counter = 0;
#label "REP" entries to go in the vcf file
if($FOREGROUND_SET eq 'REP'){
	print "Getting VDR-BVs which have discordant ASB across samples..\n";
	get_discordant_vdr_rbvs(\%discordant_vdrrbv_pos);
	my $total_discordant = scalar keys %discordant_vdrrbv_pos;
	print "Found a total of $total_discordant discordant VDR-rBV positions. VDR-rBVs at these positions will not be kept.\n";

	foreach my $var_pos (sort keys %data_fg){
		if($data_fg{$var_pos} > 1){
			$vcf_pos_counter++;
			foreach my $item (keys %fg_vcf){
				my ($chr, $pos) = (split /\t/, $item)[0,1];
				my $this_var_pos = $chr . '-' . $pos;
				if($this_var_pos eq $var_pos){
					unless($discordant_vdrrbv_pos{$this_var_pos}){						
						$fg_vcf_rep{$item} = 1;
					}
				}
			}
		}
	}
	print "Total (CONCORDANT + NOT CONCORDANT) REP variant positions: $vcf_pos_counter.\n";	
	
	#get numbers of lines and positions after discordant elements are removed.
	my $vcf_line_counter = keys %fg_vcf_rep;
	print "CONCORDANT REP variant lines: $vcf_line_counter\n";
	my %vcf_rep_conc_pos;
	foreach my $item (keys %fg_vcf_rep){
		my ($chr, $pos) = (split /\t/, $item)[0,1];
		my $coord = $chr . '-' . $pos;
		$vcf_rep_conc_pos{$coord} = 1;	
	}
	my $vcf_rep_conc_pos_counter = keys %vcf_rep_conc_pos;
	print "CONCORDANT REP variant positions: $vcf_rep_conc_pos_counter\n";	
}	

#save vcf of input variants, p-val thresholded and subset by all, rep, etc
#also save final vcf of input variants, deprived of mhc variants and bed of post-mhc variants
my %MHC_vars_count;
my %fg_var_pos_retained;

open (my $outstream_vcf,        q{>}, $output_fg_vcf) or die("Unable to open $output_fg_vcf: $!");
open (my $outstream_vcf_suppl,  q{>}, $output_fg_vcf_clean) or die("Unable to open $output_fg_vcf_clean: $!");
open (my $outstream_fg_list,    q{>}, $out_temp_fg_list) or die("Unable to open $out_temp_fg_list: $!");	
#vcf header
print $outstream_vcf       "##fileformat=VCFv4.1\n";
print $outstream_vcf       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
print $outstream_vcf_suppl "##fileformat=VCFv4.1\n";
print $outstream_vcf_suppl "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
if($FOREGROUND_SET eq 'REP'){
	foreach my $item (sort keys %fg_vcf_rep){
		print $outstream_vcf $item, "\n";
		
		my ($chr, $pos) = (split /\t/, $item)[0,1];
		my $coord = $chr . '-' . $pos;
		if( ($chr eq $MHC_chr) && ($pos ge $MHC_start) && ($pos le $MHC_end) ){
			$MHC_vars_count{$coord} = 1;
			next;
		}else{
			$fg_var_pos_retained{$coord} = 1;	
			print $outstream_vcf_suppl $item, "\n";
		}
	}	
	my $vcf_rep_conc_pos_nomhc_counter = keys %fg_var_pos_retained;
	my $MHC_count = keys %MHC_vars_count;
	print "FG: $MHC_count variant positions in the MHC discarded\n";
	print "FG: CONCORDANT,NO MHC REP variant positions retained: $vcf_rep_conc_pos_nomhc_counter\n";	
	
	foreach my $item (keys %fg_var_pos_retained){
		my ($chr, $pos) = split('-', $item);
		print $outstream_fg_list $chr . "\t" . ($pos - 1) . "\t" . $pos . "\n";		
	}
}else{
	foreach my $item (sort keys %fg_vcf){
		print $outstream_vcf $item, "\n";
		
		my ($chr, $pos) = (split /\t/, $item)[0,1];
		my $coord = $chr . '-' . $pos;
		if( ($chr eq $MHC_chr) && ($pos ge $MHC_start) && ($pos le $MHC_end) ){
			$MHC_vars_count{$coord} = 1;
			next;
		}else{
			$fg_var_pos_retained{$coord} = 1;
			print $outstream_fg_list $chr . "\t" . ($pos - 1) . "\t" . $pos . "\n";			
			print $outstream_vcf_suppl $item, "\n";
		}			
	}	
	my $vcf_pos_nomhc_counter = keys %fg_var_pos_retained;
	my $MHC_count = keys %MHC_vars_count;
	print "FG: $MHC_count variant positions in the MHC discarded\n";
	print "FG: NO MHC variant positions retained: $vcf_pos_nomhc_counter\n";	
	
	foreach my $item (keys %fg_var_pos_retained){
		my ($chr, $pos) = split('-', $item);
		print $outstream_fg_list $chr . "\t" . ($pos - 1) . "\t" . $pos . "\n";		
	}
}
close $outstream_fg_list;
close $outstream_vcf;
close $outstream_vcf_suppl;

################
#1.2 get a bed with the background variants. Remove indels
################
tie *FILE,   'IO::Zlib', $in_background, "rb";
while (<FILE>)	{ 
	next if($_ eq '');	
	next if($_ =~ /^\#/);
	chomp;
	my ($chr, $pos, $id, $info) = (split /\t/)[0,1,2,7];	
	
	#remove cases which have "esv" instead of rs id
	next if($id =~ /^esv/);
	#remove INDELS
	next if($info =~ /INDEL/);
	my $bg_var = $chr . '-' . $pos; 
	$data_bg{$bg_var} = 1; 
}
close FILE;

%MHC_vars_count = ();
my %bg_var_pos_retained;
open (my $outstream_bg_list,  q{>}, $out_temp_bg_list) or die("Unable to open $out_temp_bg_list: $!");		
foreach my $item (sort keys %data_bg){
	my ($chr, $pos) = split("-", $item);
	my $coord = $chr . '-' . $pos;
		
	if( ($chr eq $MHC_chr) && ($pos ge $MHC_start) && ($pos le $MHC_end) ){
		$MHC_vars_count{$coord} = 1;
		next;
	}else{
		$bg_var_pos_retained{$coord} = 1;
	}
}

my $vcf_pos_nomhc_counter = keys %bg_var_pos_retained;
my $MHC_count = keys %MHC_vars_count;
print "BG: $MHC_count variant positions in the MHC discarded\n";
print "BG: NO MHC variant positions retained: $vcf_pos_nomhc_counter\n";	

foreach my $item (keys %bg_var_pos_retained){
	my ($chr, $pos) = split('-', $item);
	print $outstream_bg_list $chr . "\t" . ($pos - 1) . "\t" . $pos . "\n";		
}
close $outstream_bg_list;


################
#2 create a hash to get all the fg/bg variants which are NOT in an LD block
################
system "$BEDTOOLS intersect -c -a $out_temp_fg_list -b $INPUT_LD > $out_temp_fg_buffer";
system "$BEDTOOLS intersect -c -a $out_temp_bg_list -b $INPUT_LD > $out_temp_bg_buffer";
store_variants_not_in_LDblock($out_temp_fg_buffer, \%data_fg_nold);
store_variants_not_in_LDblock($out_temp_bg_buffer, \%data_bg_nold);

###################
#3 map LD block with at least a fg/bg variant in them to variant coordinate
###################
system "$BEDTOOLS intersect -wo -a $INPUT_LD -b $out_temp_fg_list > $out_temp_fg2_buffer";
system "$BEDTOOLS intersect -wo -a $INPUT_LD -b $out_temp_bg_list > $out_temp_bg2_buffer";
map_LDblocks_to_variants($out_temp_fg2_buffer, \%data_ld2fgvar);
map_LDblocks_to_variants($out_temp_bg2_buffer, \%data_ld2bgvar);

###################
#4 [FOREGROUND ONLY] get per-sample affinity binding phenotype info (to find the best VDR-hcBV in an LD block)
###################
get_binding_affinity_levels(\%position2sample2readdepth);

####################
#5 go through all the LD blocks intersecting FOREGROUND variables and choose 1 foreground var for each
####################
foreach my $ld_block (keys %data_ld2fgvar){
	my @fgvars_in_ld; 
	my $candidate_fgvar;
	
	foreach my $fg_var_in_ld (sort keys %{ $data_ld2fgvar{$ld_block} }  ){
		push(@fgvars_in_ld, $fg_var_in_ld);
	}

	if(@fgvars_in_ld == 1){
		$candidate_fgvar = $fgvars_in_ld[0];
	}else{
		#call routine to select the "best" and fill candidate_vdrucbv
		$candidate_fgvar = get_best_fg_var(@fgvars_in_ld);
	}

	if($candidate_fgvar){
		$candidate_foreground_vars{$candidate_fgvar} = 1;
	}else{
		print STDERR "WARNING: candidate foreground variant undefined. Skipping..\n";
		next;
	}
}



####################
#6 go through all the LD blocks intersecting BACKGROUND variants and choose  one for each
#if $BACKGROUND_SET is 'R' choose a random one
#if instead $BACKGROUND_SET is 'NR' choose the foreground one in the same LD block if available, otherwise a random one
####################
foreach my $ld_block (keys %data_ld2bgvar){
	my @bg_vars_in_ld; 
	my $candidate_bgvar;
	my $fg_var_match;
	
	foreach my $bg_var_in_ld (sort keys %{ $data_ld2bgvar{$ld_block} }  ){
		push(@bg_vars_in_ld, $bg_var_in_ld);
		
		if( ($BACKGROUND_SET eq 'NR') && ($candidate_foreground_vars{$bg_var_in_ld}) ){
			$fg_var_match = $bg_var_in_ld;
		}
	}
	
	if(@bg_vars_in_ld == 1){
		$candidate_bgvar = $bg_vars_in_ld[0];
	}else{ #more than one variant in the block, pick one.
		if($BACKGROUND_SET eq 'R'){ #pick at random
			$candidate_bgvar = $bg_vars_in_ld[rand @bg_vars_in_ld];
		}elsif($BACKGROUND_SET eq 'NR'){ #pick the foreground one if available, otherwise pick at random
			if($fg_var_match){
				$candidate_bgvar = $fg_var_match;
			}else{
				$candidate_bgvar = $bg_vars_in_ld[rand @bg_vars_in_ld];
			}
		}else{
			print STDERR "Error: background set flag: $BACKGROUND_SET not recognised. Aborting..\n";
			exit -1;
		}

	}

	if($candidate_bgvar){
		$candidate_background_vars{$candidate_bgvar} = 1;
	}else{
		print STDERR "WARNING: candidate BACKGROUND variable undefined. Skipping..\n";
		next;
	}
}


###############
#7 save the foreground and background variables which are not in LD blocks
###############
foreach my $item (sort keys %data_fg_nold){
	$candidate_foreground_vars{$item} = 1;
}
foreach my $item (sort keys %data_bg_nold){
	$candidate_background_vars{$item} = 1;
}
#############
#8 print foreground and background
#############
my $tot_fg_cand = scalar keys %candidate_foreground_vars;
my $tot_bg_cand = scalar keys %candidate_background_vars;
print "After LD-clumping, $tot_fg_cand FG variant positions are retained.\n";
print "After LD-clumping, $tot_bg_cand BG variant positions are retained.\n";

print_output_bed($output_fg, \%candidate_foreground_vars);
print_output_bed($output_bg, \%candidate_background_vars);
#get vcf of the output (you need it in do_funseq_ET_02_freqmatch_bootstrap.pl)
print "Creating background vcf.gz file..\n";
system "$BEDTOOLS intersect -a $INPUT_VARIANTS -b $output_bg | sort -k1,1V -k2,2g | uniq | gzip -c > $output_bg_vcf";
print "Cleaning up..\n";
unlink $out_temp_bg_buffer; unlink $out_temp_bg2_buffer;
unlink $out_temp_fg_buffer; unlink $out_temp_fg2_buffer;
print "FINISHED.\n";























#-----------------------------------------------------------------------------------------------------------
#subroutines
#-----------------------------------------------------------------------------------------------------------
sub print_output_bed{
	my ($file, $hash) = @_;	
	
	open (my $outstream,  q{>}, $file) or die("Unable to open $file: $!");	
	foreach my $item (sort keys %$hash){
		my ($chr, $pos) = split('-',$item);
		print  $outstream $chr . "\t" . ($pos - 1) . "\t" . $pos . "\n";
	}
	close $outstream;	
	return;
}


sub store_variants_not_in_LDblock{
	my ($file, $hash) = @_;
	
	#format:
	#22	27068481	27068482	1
	#22	28807624	28807625	1
	#22	37258266	37258267	1
	#22	37258266	37258267	0  <- I want a hash with these
	open (my $instream,  q{<}, $file) or die("Unable to open $file: $!");		
	while(<$instream>){
		chomp;
		next if($_ eq '');	
		my ($chr, $pos, $in_ld) = (split /\t/)[0,2,3];
		if($in_ld == 0){
			my $coord = $chr . '-' . $pos;
			$$hash{$coord} = 1;
			next;
		}
	}
	close $instream;
	return;	
}


sub map_LDblocks_to_variants{
	my ($file, $hash) = @_;
	
	#format:
	#chr	ld_start	ld_end	info_snps	vdrhcbv_chr	vdrhcbv_start	vdrhcbv_stop	n_bases
	#1       43230559        43232798        (2)rs2816594|rs11210712 2.24    1       43232797        43232798        1
	open (my $instream,  q{<}, $file) or die("Unable to open $file: $!");		
	while(<$instream>){
		chomp;
		next if($_ eq '');	
		my ($ld_chr, $ld_start, $ld_stop, $var_chr, $var_stop) = (split /\t/)[0,1,2,5,7];
		unless($ld_chr eq $var_chr){
			print STDERR "map_LDblocks_to_variants() - Error: chromosomes do not coincide: $ld_chr vs $var_chr. Aborting..\n";
			return undef;
		}
		my $chr =  $ld_chr;
		my $ld_coord = $chr . '-' . $ld_start . '-' . $ld_stop;
		my $var_coord = $chr . '-' . $var_stop;
		$$hash{$ld_coord}{$var_coord} = 1;
	}
	close $instream;
	return;
}

#chr-position => samplename => alleles / read counts info
sub get_binding_affinity_levels{
	my ($hash) = @_;
	
	open (my $instream,  q{<}, $INPUT_ALLELESEQ_RAW) or die("Unable to open $INPUT_ALLELESEQ_RAW : $!");
	
	while(<$instream>){
		chomp;
		next if($_ =~ /^sample/); #header
		next if($_ eq '');
		
		#format
		#sample	chrm	snppos	ref	mat_gtyp	pat_gtyp	c_gtyp	phase	mat_all	pat_all	cA	cC	cG	cT	winning	SymCls	SymPval	BindingSite	cnv
		#NA06986	1	1080920	G	S	W	R	PHASED	G	A	2	0	7	0	M	Sym	0.1796875	1	1.0
		my ($sample_id, $chr,$pos, $ref, $m_allele, $p_allele, $cA, $cC, $cG, $cT, $winner, $symcls) = (split /\t/)[0,1,2,3,8,9,10,11,12,13,14,15];
		
		next if(!$chr);
		next if(!$pos);
		next unless ($symcls =~ /Asym/);
		#$chr = 'chr' . $chr;
		my $coordinate_id =  $chr . '-' . $pos; #b37
		
		my $read_data = join(",", $ref, $cA, $cC, $cG, $cT);
		$$hash{$coordinate_id}{$sample_id} = $read_data;
	}
	close $instream;
	return;
}

#sample	chrm	snppos	ref	mat_gtyp	pat_gtyp	c_gtyp	phase	mat_all	pat_all	cA	cC	cG	cT	winning	SymCls	SymPval	BindingSite	cnv
#good
#NA19190	14	88472612	C	W	M	Y	PHASED	T	C	0	6	0	0	P	Asym	0.03125	0	1.0
#NA19236	14	88472612	C	W	Y	Y	PHASED	T	C	0	7	0	0	P	Asym	0.015625	0	1.0
#not good
#NA12872	16	28857645	A	G	A	R	PHASED	G	A	9	0	22	0	M	Asym	0.0294493734837	1	1.0
#NA19190	16	28857645	A	M	G	R	PHASED	A	G	6	0	0	0	M	Asym	0.03125	1	1.0
#also not good
#NA10847	16	11766905	A	W	G	R	PHASED	A	G	0	0	6	0	P	Asym	0.03125	0	1.0
#NA19189	16	11766905	A	W	G	R	PHASED	A	G	6	0	0	0	M	Asym	0.03125	0	1.0
#the routine below will save in a hash cases where, for the same position, 'winning' is different across samples or, if winning is the same, check that the best covered allele is the same
sub get_discordant_vdr_rbvs{
	my ($hash) = @_;
	my %temp_hash;

	open (my $instream,  q{<}, $INPUT_ALLELESEQ_RAW) or die("Unable to open $INPUT_ALLELESEQ_RAW : $!");
	while(<$instream>){
        	chomp;
                next if($_ =~ /^sample/); #header
                next if($_ eq '');

		my ($sample_id, $chr,$pos, $cA, $cC, $cG, $cT, $winner, $symcls) = (split /\t/)[0,1,2,10,11,12,13,14,15];
                next if(!$chr);
                next if(!$pos);
                next unless ($symcls =~ /Asym/);

		my $coord = $chr . '-' . $pos;
		$temp_hash{$coord}{$sample_id} = join('-',$cA, $cC, $cG, $cT, $winner);
	}
	close $instream;

	#go through hash and do the checks
	foreach my $coordinate (keys %temp_hash){
		#my %winner_hash; #this must contain only one element
		my %highest_coverage_allele; #this must contain only one element
		
		foreach my $sample (keys %{ $temp_hash{$coordinate} }){
			my ($cA, $cC, $cG, $cT, $winner)  = split('-',$temp_hash{$coordinate}{$sample});
			#$winner_hash{$winner}++;
			#get highest coverage allele
			if($cA > $cC && $cA > $cG && $cA > $cT){
				$highest_coverage_allele{A}++;	
			}elsif($cC > $cA && $cC > $cG && $cC > $cT){
				$highest_coverage_allele{C}++;
			}elsif($cG > $cA && $cG > $cC && $cG > $cT){
				$highest_coverage_allele{G}++;
			}elsif($cT > $cA && $cT > $cC && $cT > $cG){
				$highest_coverage_allele{T}++;
			}else{
				print STDERR "get_discordant_vdr_rbvs(): error: two alleles have the same coverage for an ASYM vdrbv: $cA,$cC,$cG,$cT. Aborting..\n";
				exit -1;
			}
		}
		#if discordance, fill main hash
		#if(scalar keys %winner_hash > 1){
		#	$$hash{$coordinate} = 1;
		#	next;
		#}
		#there could still be discordance if the winning is the same. Check actual coverage hash
		if(scalar keys %highest_coverage_allele > 1){
			#TODO for now just discard. If you lose too many, check the proportion per allele
			$$hash{$coordinate} = 1;
            next;
		}
	}
	return;
}


sub get_best_fg_var{
	my (@ld_variant_vector) = @_;
	my %vdrbv_to_magnitude_hash;
	
	#need to cycle by sample
	foreach my $vdrbv (@ld_variant_vector){
		my @persample_dba_magns;
		#print "\n\nVDRBV: " . $vdrbv . "\n";
		
		#check if VDR-BV is VDR-QTL (it might not be in the ASB hash)
		unless(defined $position2sample2readdepth{$vdrbv}){
			#then this is a QTL. There should be only a dozen or so in the VDR-hcBV set. Keep it?
			print STDERR "WARNING: get_vdrucbv(): VDR-hcBV: $vdrbv is not in the VDR-ASB set so it's a VDR-QTL.\n";
			print STDERR "Coverage unknown: assigning magnitude 0.\n";
			$vdrbv_to_magnitude_hash{$vdrbv} = 0;
			next;
		}
		
		foreach my $sample (sort keys %{ $position2sample2readdepth{$vdrbv} }  ){
				my ($ref, $cA,$cC,$cG,$cT) = split(",", $position2sample2readdepth{$vdrbv}{$sample});
				if( ($cA eq 0) && ($cC eq 0) && ($cG eq 0) && ($cT eq 0) ){
					print STDERR "get_vdrucbv(): error: 0 read coverage for all possible nucleotides at this event ($cA,$cC,$cG,$cT). Aborting..\n";
					exit -1;
				}
				my @allele_rd = ($cA, $cC, $cG, $cT);
				my @sorted_allele_rd = sort {$b <=> $a}  @allele_rd;
				#print 'sample: ' . $sample . ': ' . join(",",@sorted_allele_rd) . "\n";
				
				#pop the largest
				my $max_coverage_allele = shift(@sorted_allele_rd);
				if($max_coverage_allele == 0){
					print STDERR "WARNING: get_vdrucbv(): max coverage allele is 0 in this sample. Skipping..\n";
					next;
				}
				#pop the second largest
				my $second_coverage_allele = shift(@sorted_allele_rd);
				
				#if a third is non zero, skip this sample altogether (mismapping)
#				my $MISMATCH_FLAG;
#				foreach my $item (@sorted_allele_rd){
#					if($item != 0){
#						print STDERR "WARNING: get_vdrucbv(): more than 2 alleles have coverage: $item. Skipping this sample..\n";
#						$MISMATCH_FLAG = 1;
#					}
#				}
#				next if($MISMATCH_FLAG);
				
				my $dba_magnitude_thissample = $max_coverage_allele - $second_coverage_allele;
				push (@persample_dba_magns, $dba_magnitude_thissample);
		}
		#print 'Vector of per-sample magnitudes: ' .  join(',', @persample_dba_magns) . "\n";
		if(@persample_dba_magns == 1){
			$vdrbv_to_magnitude_hash{$vdrbv} = $persample_dba_magns[0];
		}else{
			my $dba_magnitude = mean(@persample_dba_magns);
			#print 'Mean per-sample magnitude: ' .  $dba_magnitude, "\n";
			$vdrbv_to_magnitude_hash{$vdrbv} = $dba_magnitude;
		}
		
	}
	#return the hash key with the largest value;
	return reduce { $vdrbv_to_magnitude_hash{$a} > $vdrbv_to_magnitude_hash{$b} ? $a : $b } keys %vdrbv_to_magnitude_hash;
	#print '------------VDR-hucBV: ' . $out . "\n";
}

#------mean
sub mean {
    return sum(@_)/@_;
}
