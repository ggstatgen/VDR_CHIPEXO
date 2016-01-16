#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#23/3/2015 Chris wants me to discard all positions that do not have ancestral state. Maybe count them to see how many are there?

#I want to build a model of how the ASB-snps provoking a MOTIFBR in VDRRXR correlate with the phenotype. 
#I want to see if the GENOTYPE CHANGE and the PHENOTYPE change are in agreement

#I want to be able to restrict the VDR:RXR motif breaks to those happening at enhancers, or class II sites, or class I sites
#for the latter two cases, I need to input a bed file to intersect.

#I can do this using, for the gt, the ratio of motif break scores and for the ph, the ratio of read depths

#eg
#                    ____T_____ RC_T = 5
#                    ____C_____ RC_C = 2
#                    ____C_____
#                    ____T_____
#                    ____T_____
#                    ____T_____
#                    ____T_____
#---------------------------------------------------

#.....mat sequence.....T[T]CA..........(reference) [MOTIF BREAK SCORE = MTBS_r = 10]
#.....pat sequence.....T[C]CA..........(reference) [MOTIF BREAK SCORE = MTBS_a = 1]
#gt_score = log [ MTBS_r / MTBS_a ]
#ph_score = log [ RC_T / RC_C ] 

#I WANT TO SEE IF THE FOLLOWING IS TRUE
#gt_score > 0 && ph_score > 0
#or
#gt_score < 0 && ph_score < 0 

#Cases of the kind:
#gt_score > 0 && ph_score < 0
#or
#gt_score > 0 && ph_score < 0
#ARE BAD: disagreements 

#You need to identify the ancestral reference both from the [RAW read counts] and from the [funseq motif break] (consult 1kg hash BOTH times)
#Input:
#1)######################################################
#Output.vcf file from funseq - I need lines tagged with MOTIFBR and VDR_JASPAR

#eg
#MOTIFBR=MAX#Myc_known9_8mer#102248644#102248656#-#9#0.068966#0.931034
#TF name # motif name # motif start # motif end # motif strand # mutation position # alternative allele frequency in PWM # reference allele frequency in PWM
#MOTIFG=GATA_known5#75658824#75658829#-#1#4.839#4.181
#motif name # motif start # motif end # motif strand # mutation position # sequence score with alternative allele # sequence score with reference allele. (0-based, end exclusive)

#2)######################################################
#all raw alleleseq outputs with raw count reads from REF carrying allele parent  and ALT-carrying allele parent, for each sample
#typical input
#chrm	snppos   	ref	mat_gtyp	pat_gtyp	c_gtyp	phase	mat_all	pat_all	cA	cC	cG	cT	winning	SymCls	SymPval	BindingSite	cnv
#1	121428238	C	M	Y	M	PHASED	A	C	13	3	0	0	M	Asym	0.021270752	0	1
#1	121428255	T	M	Y	W	PHASED	A	T	19	0	0	0	M	Asym	3.81469726562E-006	0	1

#read fields are
#8,9,10,11,12,13

#3)######################
#ALL REFERENCES to REF ALLELE mean ANCESTRAL ALLELE to me (given my funseq analysis)
#So the third input is the full set of 1kg snps tested for ASB/QTL, from which I will get ancestral information

#for each motifbr position, for eack sample, you will gather a PAIR (fc_genotype, fc_phenotype)
#you need to distinguish the cases [ (+,+) and (-,-) ] from [(+,-) and (-,+) ]

my $infile;
my $INFILE_ALLELESEQ_RAW = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/_RAW/interestingHets_vdrchipexo.txt";
my $INPUT_VARIANTS = '/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/ALLELESEQ_20samples_merged.vcf.gz';
my $motif_name = 'VDR_JASPAR';
my $ENH_ONLY;
my $subset_bed; #if you only want to collect data for positions in an external bed file (eg CLASS I sites)
my $sym_variants; #flag, set it if you're dealing with sym variants
GetOptions(
        'i=s'      =>\$infile,
        'subset=s' =>\$subset_bed,
        'sym'      =>\$sym_variants,
        'enh'      =>\$ENH_ONLY
); 
#$infile = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/Output.vcf";
#$subset_bed = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_motif_scan_pscanchip/Pscanchip_occurrences_VDRRXR_jaspar_classI.bed";

if(!$infile){
	print "USAGE: do_funseq_collect_MOTIFBR_gtph.pl -i=<INFILE> -subset=<SUBSET_BED> -sym -enh\n";
    print "<INFILE> vcf output of funseq2\n";
    print "(optional)<SUBSET_BED> only collect data for VDR-BVs breaking motifs intersecting with this bed file(default=no)\n";
    print "(optional)<sym> flag; set it if you are dealing with SYM variants from alleleseq (default=no)\n";
    print "(optional)<enh> whether to only use snps with NCENC=enhancer (default=no)\n";
    exit 1;
}
my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
my $output;

if($subset_bed){
	$output  = $directory . $basename . '_concordance_subsetbed.Rdata';
}elsif($ENH_ONLY){
	$output  = $directory . $basename . '_concordance_enhonly.Rdata';
}else{
	$output  = $directory . $basename . '_concordance_all.Rdata';
}

################
#0 if there is a bed file with peaks (eg class I peaks) save their coordinates in hash
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
#1 build huge hash to map chr-pos to ref,ancestral and alternate
#####################
my %variants_1kg;
tie *FILE,   'IO::Zlib', $INPUT_VARIANTS, "rb";
while (<FILE>)	{ 
	chomp;
	next if($_ eq '');
	next if($_ =~ /^#/);
	
	my @fields = split("\t", $_);
	#skip indels
	next if(length($fields[3]) > 1);
	next if(length($fields[4]) > 1);
	
	my $key = $fields[0] . '-' . $fields[1];
	my $ref = uc($fields[3]); 
	my $alt = uc($fields[4]);
	
	#get ancestral
	my @info = split(";", $fields[7]);	
	if($info[0] =~ /^AA=(.*)/){
		my $anc = uc($1);
		if($anc eq $alt){
			$variants_1kg{$key}{ANC} = $alt;
			$variants_1kg{$key}{DER} = $ref;	
		}elsif($anc eq $ref){
			$variants_1kg{$key}{ANC} = $ref;
			$variants_1kg{$key}{DER} = $alt;			
		}elsif( ($anc eq '') or ($anc eq 'N') or ($anc eq '.') or ($anc eq '-') ){
			#$variants_1kg{$key}{ANC} = 'NA';
			#$variants_1kg{$key}{DER} = 'NA';
			next;
		}else{ 
			#the ancestral allele is a base, but a different one from ref and alt
			#so both ref and alt are "derived"
			#$variants_1kg{$key}{ANC} = 'NA';
			#$variants_1kg{$key}{DER} = 'NA';
			#Question asked to 1kg staff
			#Ignore these cases. The primates all had an allele and then a variation happened
			next;		
		}
	}else{
			#$variants_1kg{$key}{ANC} = 'NA';
			#$variants_1kg{$key}{DER} = 'NA';
			next;	
	}	
	
#	if($info[0] =~ /^AA=(.*)/){
#		my $anc = uc($1);
#		if($anc eq $alt){
#			$variants_1kg{$key}{REF} = $alt;
#			$variants_1kg{$key}{ALT} = $ref;			
#			next;
#		}
#	}
#	$variants_1kg{$key}{REF} = $ref;
#	$variants_1kg{$key}{ALT} = $alt;
}
close FILE;
#This should have taken a long time. now you know what the ancestral reference and alternate reference is.

###################
#2 build interestinghets hash
###################
my %position2sample2readdepth;
open (my $instream,  q{<}, $INFILE_ALLELESEQ_RAW) or die("Unable to open $INFILE_ALLELESEQ_RAW : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^sample/); #header
	next if($_ eq '');
	
	#format
	#sample	chrm	snppos	ref	mat_gtyp	pat_gtyp	c_gtyp	phase	mat_all	pat_all	cA	cC	cG	cT	winning	SymCls	SymPval	BindingSite	cnv
	#NA06986	1	1080920	G	S	W	R	PHASED	G	A	2	0	7	0	M	Sym	0.1796875	1	1.0
	my ($sample_id, $chr,$snppos, $ref, $m_allele, $p_allele, $cA, $cC, $cG, $cT, $winner, $symcls) = (split /\t/)[0,1,2,3,8,9,10,11,12,13,14,15];
	next if(!$chr);
	next if(!$snppos);
	
	if($sym_variants){
		next unless ($symcls =~ /Sym/);
	}else{
		next unless ($symcls =~ /Asym/);
	}
	
	my $coordinate_id = $chr . '-' . $snppos;
	my $allele_ancestral = $variants_1kg{$coordinate_id}{ANC};
	my $allele_derived = $variants_1kg{$coordinate_id}{DER}; 
	next unless($allele_ancestral);
	next unless($allele_derived);
	
	#if($ref ne $ancestral_ref){
	#	print STDERR "$coordinate_id: ref and ancestral differ. ref: $ref, ancestral_ref: $ancestral_ref, ancestral_alt: $ancestral_alt\n";
	#}
	
	if( ($m_allele eq 'None') or ($p_allele eq 'None')  ){
	}else{
		if( ($allele_ancestral eq $m_allele)  && ($allele_derived eq $p_allele) ){		
		}elsif( ($allele_ancestral eq $p_allele)  && ($allele_derived eq $m_allele)  ){	
		}else{
			print STDERR "$coordinate_id: (ancestral ref/alt and mat/pat don't match: ($allele_ancestral, $allele_derived) and ($m_allele, $p_allele)\n";
			next;
		}		
	}
	#read data string format: ref,m_allele,p_allele,$cA,$cC,$cG,$cT,winner;
	#my $read_data = join(",", $ref,$m_allele, $p_allele, $cA, $cC, $cG, $cT, $winner);
	#read data string format: ancestral_ref, ancestral_alt, $cA,$cC,$cG,$cT;
	my $read_data = join(",", $allele_ancestral,$allele_derived, $cA, $cC, $cG, $cT);
	$position2sample2readdepth{$coordinate_id}{$sample_id} = $read_data;
}
close $instream;


#########################
#4 pick asb/sym snps breaking VDR:RXR  motifs from funseq
#########################
my %motifmodel_motifpos2geomicpos2concordance;
open ($instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	my $info_motifbr;
	my $info_ncenc;
	#my $info_motifg; I don't think motifgain was working 
	next if($_ =~ /^\#/);
	next if($_ eq '');
	if($ENH_ONLY){ next unless($_ =~ /NCENC/); }
	
	#get info field
	my ($chr, $pos, $rs_id, $funseq_ref, $funseq_alt, $info) = (split /\t/)[0,1,2,3,4,7];
	my $chr_hg19 = $chr;
	$chr =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19
	
	my @info = split(";", $info);
	foreach my $item (@info){
		$info_motifbr = $item if($item =~ /^MOTIFBR/);
		#$info_motifg = $item if($item =~ /^MOTIFG/);
		$info_ncenc = $item if($item =~ /^NCENC/);
	}
	next unless($info_motifbr);
	
	#EHNANCER FILTER
	if($ENH_ONLY){ next unless($info_ncenc =~ /Enhancer/); }
	#BED FILTER
	if($subset_bed){ next unless(defined check_coords_in_bed($chr_hg19, $pos)); }
	my $genomic_coord = $chr . '-' . $pos;
	
	#motifbr info is structured as follows
	#MOTIFBR=VDR#VDR_JASPAR#139340766#139340781#+#13#0.000#0.900,VDR#VDR_XXmotif#139340766#139340781#+#13#0.04338#0.80145
	my ($motif_change_type,$motif_change_byTF) = split("=", $info_motifbr);
	my @motif_change_byTF = split(",", $motif_change_byTF);
	
	foreach my $item (@motif_change_byTF){
		my ($TF,$TF_motif,$TF_motif_start,$TF_motif_end,$TF_motif_strand,$TF_snp_position,$TF_alt_score,$TF_ref_score) =  split("#", $item);		
		next unless( $TF_motif eq $motif_name);
		
		if(!$position2sample2readdepth{$genomic_coord}){
			print STDERR "$rs_id - $genomic_coord: no ancestral/derived info available for this VDR-BV. Skipping..\n";
			next;
		}		
		
		my $GT_SCORE;
		my $concordant = 0;
		my $discordant = 0;
		#my @concordance_vector;
		#Compute gt_score from $TF_alt_score and $TF_alt_score
		#these should already account for ancestral allele in funseq
		$GT_SCORE = ( $TF_ref_score + 1 ) / ($TF_alt_score + 1); 
		
		#Compute ph score and evaluate concordance with gt score at that position	
		foreach my $sample (sort keys %{ $position2sample2readdepth{$genomic_coord} }){
			my $PH_SCORE_this_sample = process_read_data($position2sample2readdepth{$genomic_coord}{$sample}, $sample, $rs_id);	
			if( ($GT_SCORE >= 1) && ($PH_SCORE_this_sample >= 1) ){
				#push(@concordance_vector, 'up-up');
				$concordant += 1;
			}elsif( ($GT_SCORE >= 1) && ($PH_SCORE_this_sample < 1)  ){
				#push(@concordance_vector, 'up-down');
				$discordant += 1;
			}elsif( ($GT_SCORE < 1) && ($PH_SCORE_this_sample >= 1)  ){
				#push(@concordance_vector, 'down-up');
				$discordant += 1;
			}elsif( ($GT_SCORE < 1) && ($PH_SCORE_this_sample < 1)  ){
				#push(@concordance_vector, 'down-down');
				$concordant += 1;
			}else{
				print STDERR "Should never be here\n";
			}
		}
		
		#attach the concordance counts to the main hash, after having adjusted the break by strand:
		my $break_position;
		my $motif_size = $TF_motif_end - $TF_motif_start;
		if($TF_motif_strand =~ /-/){
			my $TF_snp_position_RC = ($motif_size - $TF_snp_position + 1 );
			$break_position = 'pos_' . $TF_snp_position_RC;
		}elsif($TF_motif_strand =~ /\+/){
			$break_position = 'pos_' . $TF_snp_position;
		}else{
			print STDERR "MOTIFBREAK strand: $TF_motif_strand not recognised, skipping..\n";
			next;
		}
		#save the concordance data for this motif position at this genomic coordinate
		$motifmodel_motifpos2geomicpos2concordance{$break_position}{$genomic_coord}{AGREE} = $concordant;
		$motifmodel_motifpos2geomicpos2concordance{$break_position}{$genomic_coord}{DISAGREE} = $discordant;
	}
}
close $instream;

open (my $outstream,  q{>}, $output) or die("Unable to open $output : $!");
print $outstream "MOTIF_POS\tCONCORDANT\tDISCORDANT\n";
foreach my $motif_pos (sort keys %motifmodel_motifpos2geomicpos2concordance){
	my $sum_concordant_for_this_motif_pos = 0;
	my $sum_discordant_for_this_motif_pos = 0;
	
	foreach my $genomic_pos (sort keys %{ $motifmodel_motifpos2geomicpos2concordance{$motif_pos} }){
		$sum_concordant_for_this_motif_pos += $motifmodel_motifpos2geomicpos2concordance{$motif_pos}{$genomic_pos}{AGREE};
		$sum_discordant_for_this_motif_pos += $motifmodel_motifpos2geomicpos2concordance{$motif_pos}{$genomic_pos}{DISAGREE};
	}
	print $outstream "$motif_pos\t$sum_concordant_for_this_motif_pos\t$sum_discordant_for_this_motif_pos\n";
}
close $outstream;

#subroutines===============================================================================

#1
#this splits a read_data string of the kind "c,a,c,0,6,0,0,P" or "g,none,none,0,0,0,6"
#into a phenotype disruption score lnFC = ln [ (anc_count + 1) / (alt_count + 1) ] (or lnFC = ln [ (ref_count + 1) / (alt_count + 1) ] )
#if the motif break snps has lowered the coverage at the position, this subroutine should output a positive value

#input meaning:
##read data string format: ref,m_allele,p_allele,$cA,$cC,$cG,$cT,winner;
#read data string format: ancestral_ref, ancestral_alt, $cA,$cC,$cG,$cT;

sub process_read_data{
	my ($read_data_string, $lcl_sample, $id) = @_;
	
	my ($allele_ancestral,$allele_derived,$cA,$cC,$cG,$cT) = split(",", $read_data_string);
	my $anc_count; my $der_count;
	my $anc_base = $allele_ancestral;
	my $der_base = $allele_derived;

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
	my $fold_change = ( $anc_count + 1 ) / ( $der_count + 1 );
	return $fold_change;
}

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
