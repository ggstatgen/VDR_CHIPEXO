#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;


#30.5.2016
#This script is an extension of do_REVIEW_get_bestscoring_PWM_at_VDRBVs.pl that looks at concordance before assigning VDR-BVs to motifs.

#22/5/2016
#Following discussion with Chris (few days before leaving Oxford) he reckons that if we extend this analysis to clean up based on phenotype
# (concordance and discordance) we might be able to create a response for the reviewers.
#Therefore I want to extend this script to only consider VDR-BVs under VDR:RXR which are concordant in impact on genotype and impact on phenotype
#to calculate concordance, we will need:
#1-genotypes, and ancestral info, so the full alleleseq 20 sample vcf
#2-phenotypes and therefore the alleleseq  output

#see  do_REVIEW_get_bestscoring_PWM_at_VDRBVs.pl for full readme


my $BEDTOOLS = `which bedtools`; chomp $BEDTOOLS;
my $RSCRIPT = `which RRscript`; chomp $RSCRIPT;
my $IN_VDRBV = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_REVIEW_asym_anc/Output_noDBRECUR.vcf";
#my $IN_VDRBV = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/SUPPL_DATA_Output_noDBRECUR_REP_hg19.vcf";
my $PWM_FILE = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Processed_PFMs_jaspar_FUNSEQ_INPUT.txt";
my $PLOT_EXT = 'pdf';

#phenotype from here
my $INFILE_ALLELESEQ_RAW = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/_RAW/interestingHets_vdrchipexo.txt";
#genotype from here
my $INPUT_VARIANTS = '/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/ALLELESEQ_20samples_merged.vcf.gz';

my $INPUT_RIS_DIR;
my $MIN_SCORE;
#GLOBALS
#1 vdr-bvs - this will be checked against the ris intervals and intersecting items will be removed
my %VDRBV_coords;
#2 contains the jaspar representation of the pws to get the correct length
my %JASPAR_MOTIF; 
#3 result structure: coord->associated_pwm->score
my %RESULTS; my %results_allpwms;
#If the best has same score for many motifs? Take the LONGEST

GetOptions(
        'i=s'		=>\$INPUT_RIS_DIR,          
        's=f'		=>\$MIN_SCORE,   
);
#temp
#$INPUT_RIS_DIR = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/VDR-BV";
#$MIN_SCORE = 0.7;

my $USAGE = "\nUSAGE: $0 -i=<VDRBV_RIS_DIR> -s=<MIN_SCORE>\n" .
			"<VDRBV_RIS_DIR> full path with .ris files from PscanChip\n" .
			"(opt)<MIN_SCORE> min PscanChip score to consider (default=undef)\n";
			
unless($INPUT_RIS_DIR){
	print $USAGE;
	exit -1;
}
print "Minimum PScanChIP score set to $MIN_SCORE\n" if ($MIN_SCORE);

#my $RXR_VDR_RIS = "Pscanchip_hg19_bkgGM12865_Jaspar_VDRBVs_RXRA-VDR_MA0074.1_sites.ris";
my $RXR_VDR_RIS = "Pscanchip_hg19_bkgGM12865_Jaspar_VDRrBVs_RXRA-VDR_MA0074.1_sites.ris";
my $RXR_VDR_PATH = $INPUT_RIS_DIR . '/' . $RXR_VDR_RIS;


#COLLECT VDR-BV genotypes
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
			next;
		}else{ 
			next;		
		}
	}else{
			next;	
	}	
}
close FILE;


#COLLECT VDR::RXR PHENOTYPES
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
			next;
		}		
	}
	my $read_data = join(",", $allele_ancestral,$allele_derived, $cA, $cC, $cG, $cT);
	$position2sample2readdepth{$coordinate_id}{$sample_id} = $read_data;
}
close $instream;

################
#3. Get the motif length from the encode representations of the motif---------------------
################
get_motif_lengths($PWM_FILE, \%JASPAR_MOTIF);

##################
#4 Build pwm ID
##################
my ($motif_string, $full_motif_id, $motif_length) = get_pwm_id($RXR_VDR_RIS);
print STDERR "The length of the motif: $full_motif_id according to the JASPAR PWM is $motif_length\n";

#########
#5. Find the motif intervals that intersect an RXR:VDR motif -> pick the *concordant* ones and set them aside.
#if the vdrbv hits a motif and the hit is concordant, save it into a result hash, and REMOVE it from the master VCF file -
#or in other words create a vcf hash with all the ones that:
#1-do not hit a VDR:RXR motifbr
#2-they do, but the phenotype is discordant
#########
my %vdrbv_breaking_rxrvdr; #these break rxrvdr with concordance
my %vdrbv_notbreaking_rxrvdr; #these break one or more motif but not rxrvdr so get tested with all other pwms
my %vdrbv_nobreak; #these do not break anything, any motif at all - count

my %motifmodel_motifpos2geomicpos2concordance;
open ($instream,  q{<}, $IN_VDRBV) or die("Unable to open $IN_VDRBV : $!");
while(<$instream>){
	chomp;
	my $info_motifbr;
	next if($_ =~ /^\#/);
	next if($_ eq '');
	
	#get info field
	my ($chr, $pos, $rs_id, $funseq_ref, $funseq_alt, $info) = (split /\t/)[0,1,2,3,4,7];
	my $chr_hg19 = $chr;
	$chr =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19
	my $genomic_coord = $chr . '-' . $pos;
	
	my @info = split(";", $info);
	foreach my $item (@info){
		$info_motifbr = $item if($item =~ /^MOTIFBR/);
	}

	#if the VDR-BV does not break any motif at all, including the enriched, keep it to count 
	unless($info_motifbr){
		$vdrbv_nobreak{$genomic_coord} = 1;
		next;
	}
	
	#motifbr info is structured as follows
	#MOTIFBR=VDR#VDR_JASPAR#139340766#139340781#+#13#0.000#0.900,VDR#VDR_XXmotif#139340766#139340781#+#13#0.04338#0.80145
	my ($motif_change_type,$motif_change_byTF) = split("=", $info_motifbr);
	my @motif_change_byTF = split(",", $motif_change_byTF);
	
	#here you only want the motif corresponding to full_motif_id
	foreach my $item (@motif_change_byTF){
		my ($TF,$TF_motif,$TF_motif_start,$TF_motif_end,$TF_motif_strand,$TF_snp_position,$TF_alt_score,$TF_ref_score) =  split("#", $item);		
		#the vdrbv breaks a motif but not the RXR:VDR
		unless( $TF_motif eq $full_motif_id){
			$vdrbv_notbreaking_rxrvdr{$_} = 1;
			next;
		}
		
		if(!$position2sample2readdepth{$genomic_coord}){
			print STDERR "$rs_id - $genomic_coord: no ancestral/derived info available for this VDR-BV. Skipping..\n";
			next;
		}		
		my $GT_SCORE;my $concordant = 0;my $discordant = 0;
		#Compute gt_score from $TF_alt_score and $TF_alt_score
		#these should already account for ancestral allele in funseq
		$GT_SCORE = ( $TF_ref_score + 1 ) / ($TF_alt_score + 1); 
		
		#Compute ph score and evaluate concordance with gt score at that position	
		foreach my $sample (sort keys %{ $position2sample2readdepth{$genomic_coord} }){
			my $PH_SCORE_this_sample = process_read_data($position2sample2readdepth{$genomic_coord}{$sample}, $sample, $rs_id);	
			if( ($GT_SCORE >= 1) && ($PH_SCORE_this_sample >= 1) ){
				$concordant += 1;
			}elsif( ($GT_SCORE >= 1) && ($PH_SCORE_this_sample < 1)  ){
				$discordant += 1;
			}elsif( ($GT_SCORE < 1) && ($PH_SCORE_this_sample >= 1)  ){
				$discordant += 1;
			}elsif( ($GT_SCORE < 1) && ($PH_SCORE_this_sample < 1)  ){
				$concordant += 1;
			}else{
				print STDERR "Should never be here\n";
			}
		}
		
		#if the VDR-BV gt change agrees with the ph, keep it, otherwise, test this vdr bv with another PWM
		if($concordant > $discordant){ #TODO rather permissive definition
			#I want to count unique coordinates so I only pass those to the hash:
			$vdrbv_breaking_rxrvdr{$genomic_coord} = 1;
		}else{
			$vdrbv_notbreaking_rxrvdr{$_} = 1;
		}
	}
}
close $instream;

my $vdrbv_nobreak = scalar keys %vdrbv_nobreak;
my $vdrbv_breaking_rxrvdr = scalar keys %vdrbv_breaking_rxrvdr;

print STDERR "Number of VDR-BVs that break none of the 30ish enriched: $vdrbv_nobreak\n";
print STDERR "Number of VDR-BVs that break RXR::VDR and are gt-ph concordant: $vdrbv_breaking_rxrvdr\n";


#now test the remaining VDR-BVs in vdrbv to identify the best
#if there are several PWMs broken, take the one with the best score (how?)
#if many PWMs with the same score are broken, take the longest
###############
#ALL other PWMs
###############
#Here I will need an intermediate hash, with
#pos -> pmw_id1 -> score1
#   |
#   |_> pwm_id2 -> score2
#For each of these positions, I will have to choose the best scoring pwm_id. Also, if two or more have the same score, I pick the largest PWM
#1st pass: label each vdrbv with the best instanes of all the pwm(s) they fall in
chdir $INPUT_RIS_DIR;
my @files = <Pscanchip_hg19*.ris>;
foreach my $RIS_FILE (@files){
	next if ($RIS_FILE eq $RXR_VDR_RIS);
	my %concordant_for_this_pwm;
	
	my $RIS_FILE_PATH = $INPUT_RIS_DIR  . '/' . $RIS_FILE;
	
	#1. Build pwm ID------------------------------------------------------------------------------
	my ($motif_string, $full_motif_id, $motif_length) = get_pwm_id($RIS_FILE);
	#print STDERR "The length of the motif: $full_motif_id according to the JASPAR PWM is $motif_length\n";
	
	#2. Write ris into bed------------------------------------------------------------------------
	my $tmp_ris_bed           = $INPUT_RIS_DIR . '/TMP_from_ris_'     . $motif_string . '.bed';
	my $tmp_vdrbv_concord_bed = $INPUT_RIS_DIR . '/TMP_vdrbvconc_'    . $motif_string . '.bed';
	my $tmp_intersect_bed     = $INPUT_RIS_DIR . '/TMP_intersect_'    . $motif_string . '.bed';
	write_ris_to_bed_file($RIS_FILE_PATH, $tmp_ris_bed, $motif_length);
	
	#3. Test all the VDR-BVs not breaking RXR::VDR for concordant breakage of this full motif id---------------
	foreach (keys %vdrbv_notbreaking_rxrvdr){
		my $info_motifbr;
		my ($chr, $pos, $rs_id, $funseq_ref, $funseq_alt, $info) = (split /\t/)[0,1,2,3,4,7];
		my $chr_hg19 = $chr;
		$chr =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19
		my $genomic_coord = $chr . '-' . $pos;
		my $bedline = $chr . "\t" . ($pos-1) . "\t" . $pos;
		
		my @info = split(";", $info);
		foreach my $item (@info){
			$info_motifbr = $item if($item =~ /^MOTIFBR/);
		}
		
		#double check - it should never enter the following 
		unless($info_motifbr){
			print STDERR "WARNING: it shouldn't be here...Skipping this vdrbv..\n";
			next;
		}
		
		#motifbr info is structured as follows
		#MOTIFBR=VDR#VDR_JASPAR#139340766#139340781#+#13#0.000#0.900,VDR#VDR_XXmotif#139340766#139340781#+#13#0.04338#0.80145
		my ($motif_change_type,$motif_change_byTF) = split("=", $info_motifbr);
		my @motif_change_byTF = split(",", $motif_change_byTF);
		
		#here you only want the motif corresponding to full_motif_id
		foreach my $item (@motif_change_byTF){
			my ($TF,$TF_motif,$TF_motif_start,$TF_motif_end,$TF_motif_strand,$TF_snp_position,$TF_alt_score,$TF_ref_score) =  split("#", $item);		
			#the vdrbv breaks a motif but not the one we are probing right now - skip it:
			next unless( $TF_motif eq $full_motif_id);
			
			if(!$position2sample2readdepth{$genomic_coord}){
				print STDERR "$rs_id - $genomic_coord: no ancestral/derived info available for this VDR-BV. Skipping..\n";
				next;
			}		
			
			my $GT_SCORE;my $concordant = 0;my $discordant = 0;
			#Compute gt_score from $TF_alt_score and $TF_alt_score
			#these should already account for ancestral allele in funseq
			$GT_SCORE = ( $TF_ref_score + 1 ) / ($TF_alt_score + 1); 
			
			#Compute ph score and evaluate concordance with gt score at that position	
			foreach my $sample (sort keys %{ $position2sample2readdepth{$genomic_coord} }){
				my $PH_SCORE_this_sample = process_read_data($position2sample2readdepth{$genomic_coord}{$sample}, $sample, $rs_id);	
				if( ($GT_SCORE >= 1) && ($PH_SCORE_this_sample >= 1) ){
					$concordant += 1;
				}elsif( ($GT_SCORE >= 1) && ($PH_SCORE_this_sample < 1)  ){
					$discordant += 1;
				}elsif( ($GT_SCORE < 1) && ($PH_SCORE_this_sample >= 1)  ){
					$discordant += 1;
				}elsif( ($GT_SCORE < 1) && ($PH_SCORE_this_sample < 1)  ){
					$concordant += 1;
				}else{
					print STDERR "Should never be here\n";
				}
			}
			
			#I have those that are concordant for this PWM but not the scores..I need to intersect a file of concordant vdrbvs with the pscanchip bed obtained in 2..
			$concordant_for_this_pwm{$bedline} = 1 if($concordant > $discordant);
		}
	}
	#create bed of concordant vdrbv and do the intersection down here to get the score
	open (my $outstream,  q{>}, $tmp_vdrbv_concord_bed) or die("Unable to open $tmp_vdrbv_concord_bed : $!");
	foreach my $item (keys %concordant_for_this_pwm){ print $outstream $item, "\n"; }
	close $outstream;
		
	#4. Get the pscanchip scores for the concordant vdrbvs for this pwm:-------------------------------------------------------
	print "Working on $full_motif_id:\n";
	system "$BEDTOOLS intersect -wo -a $tmp_vdrbv_concord_bed -b $tmp_ris_bed | awk -F \"\t\" \'{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$7\"\t\"\$8}' > $tmp_intersect_bed";
	#4 all of them should have intersections (TODO check) - save them with their scores
	open (my $instream,  q{<}, $tmp_intersect_bed) or die("Unable to open $tmp_intersect_bed : $!");
	while(<$instream>){
		chomp;
		next if($_ eq '');
		#get coords
		my ($chr, $pos, $score) = (split /\t/)[0,2,4];
		unless($chr =~ /^chr/){ $chr = 'chr' . $chr;}
		my $vdr_bv_coord = $chr . '-' . $pos;
		$results_allpwms{$vdr_bv_coord}{$motif_string}{'SCORE'}      = $score;
		$results_allpwms{$vdr_bv_coord}{$motif_string}{'PWM_LENGTH'} = $motif_length;
	}
	close $instream;
	unlink $tmp_ris_bed;
	unlink $tmp_intersect_bed;
	unlink $tmp_vdrbv_concord_bed;
}
#TODO WORK HERE




##2 write bed of ris intervals
#my $tmp_ris_bed               = $INPUT_RIS_DIR . '/TMP_from_ris_'    . $motif_string . '.bed';
#my $tmp_intersect_bed         = $INPUT_RIS_DIR . '/TMP_intersect_'   . $motif_string . '.bed';
#my $VDRBV_no_RXRVDR_intersect = $INPUT_RIS_DIR . '/TMP_nointersect_' . $motif_string . '.bed';
#my $tmp_intersect_vcf         = $INPUT_RIS_DIR . '/TMP_intersect_'   . $motif_string . '.vcf';
#write_ris_to_bed_file($RXR_VDR_PATH, $tmp_ris_bed, $motif_length);
#system "$BEDTOOLS intersect -wo -a $IN_VDRBV -b $tmp_ris_bed | awk -F \"\t\" \'{print \$1\"\t\"\$2-1\"\t\"\$2\"\t\"\$12\"\t\"\$13}' > $tmp_intersect_bed";
#system "$BEDTOOLS intersect -v -a $IN_VDRBV -b $tmp_ris_bed | awk -F \"\t\" \'{print \$1\"\t\"\$2-1\"\t\"\$2}\'> $VDRBV_no_RXRVDR_intersect";
#
#open ($instream,  q{<}, $tmp_intersect_bed) or die("Unable to open $tmp_intersect_bed : $!");
#while(<$instream>){
#	chomp;
#	next if($_ eq '');
#	
#	#get coords
#	my ($chr, $pos, $score) = (split /\t/)[0,2,4];
#	unless($chr =~ /^chr/){$chr = 'chr' . $chr;}
#	my $coord = $chr . '-' . $pos;
#	$RESULTS{$coord}{$motif_string} = $score;
#}
#close $instream;
#unlink $tmp_ris_bed;
#unlink $tmp_intersect_bed;











	
#2nd pass--
#for every position, choose one pwm
#criteria:
#if 1 pwm, keep it
#if 2 or more, get highest score
#if score same get longest

foreach my $vdrbv_pos (keys %results_allpwms){
	my $number_of_pwms = keys %{$results_allpwms{$vdrbv_pos}};
	
	if($number_of_pwms == 1){
		foreach my $this_pwm (keys %{$results_allpwms{$vdrbv_pos}}){
			$RESULTS{$vdrbv_pos}{$this_pwm} = $results_allpwms{$vdrbv_pos}{$this_pwm}{'SCORE'};
		}
		next;
	}
	
	my $candidate_score = 0;
	my $candidate_length = 0;
	my $candidate_pwm = '';
	foreach my $this_pwm (keys %{$results_allpwms{$vdrbv_pos}}){
		my $this_score = $results_allpwms{$vdrbv_pos}{$this_pwm}{'SCORE'};
		my $this_length =  $results_allpwms{$vdrbv_pos}{$this_pwm}{'PWM_LENGTH'};
		if($this_score > $candidate_score){
			$candidate_pwm = $this_pwm;
			$candidate_score = $this_score;
			$candidate_length = $this_length;
		}elsif($this_score == $candidate_score){
			if($this_length >= $candidate_length){
				$candidate_pwm = $this_pwm;
				$candidate_score = $this_score;
				$candidate_length = $this_length;					
			}else{
				next;
			}
		}else{
			next;
		}
	}
	$RESULTS{$vdrbv_pos}{$candidate_pwm} = $candidate_score;
}


$MIN_SCORE = 'NA' unless($MIN_SCORE);
my $out_Rdata    = $INPUT_RIS_DIR . '/RESULTS_bestscoring_pwm_scorethrs_' . $MIN_SCORE . '.Rdata';
my $out_Rcode    = $INPUT_RIS_DIR . '/RESULTS_bestscoring_pwm_scorethrs_' . $MIN_SCORE . '.R';
my $out_Rplot    = $INPUT_RIS_DIR . '/RESULTS_bestscoring_pwm_scorethrs_' . $MIN_SCORE . '.' . $PLOT_EXT;
my $output_file  = $INPUT_RIS_DIR . '/RESULTS_bestscoring_pwm_scorethrs_' . $MIN_SCORE . '.tsv';

open (my $outstream,  q{>}, $output_file) or die("Unable to open $output_file : $!");
print $outstream "CHR\tPOS\tASSIGNED_PWM\tSCORE\n";
foreach my $vdrbv_coords (keys %RESULTS){
	my ($chr, $pos) = split('-',$vdrbv_coords);
	foreach my $pwm_id (keys %{$RESULTS{$vdrbv_coords}}){
		print $outstream $chr, "\t", $pos, "\t", $pwm_id, "\t", $RESULTS{$vdrbv_coords}{$pwm_id}, "\n";
	}
}
close $outstream;
#unlink $VDRBV_no_RXRVDR_intersect;

#output proportions, too.


#R plot
#I want a bar plot, with the numbers on top of each bar
#Write R bar plots
system "grep -v \"ASSIGNED_PWM\" $output_file | cut -f 3 | sort > $out_Rdata";
open ($outstream,  q{>}, $out_Rcode) or die("Unable to open $out_Rcode : $!");
print $outstream "data <- read.table(\"$out_Rdata\",sep=\"\\t\")" . "\n";
print $outstream "data_c <- table(data)" . "\n";
print $outstream "$PLOT_EXT(file=\"$out_Rplot\", width=14,height=7)" . "\n";
print $outstream "par(mar = c(11,6,6,4) + 0.3)" . "\n";
print $outstream "bplt <- barplot(data_c, ylab=\"\#VDR-BVs\", width=1, main=\"\#VDR-BVs by best PWM model hit(t=$MIN_SCORE)\", las=2)" . "\n";
print $outstream "text(x=bplt, y=data_c, labels=as.character(data_c), pos=3, cex = 0.8, col = \"red\", xpd=TRUE)" . "\n";
print $outstream "dev.off()" . "\n";
close $outstream;

system "$RSCRIPT $out_Rcode";

















#subs--------------------------------------------------------------------------------------------
sub get_vdrbvs_from_file{
	my ($file, $hash) = @_;
	
	open (my $instream,  q{<}, $file) or die("Unable to open $file : $!");
	while(<$instream>){
		chomp;
		next if($_ =~ /^\#/);
		next if($_ eq '');
		
		#get coords
		my ($chr, $pos) = (split /\t/)[0,1];
		unless($chr =~ /^chr/){
			$chr = 'chr' . $chr;
		}
		my $coord = $chr . '-' . $pos;
		$$hash{$coord} = 1;
	}
	close $instream;	
	return;	
}


sub get_motif_lengths{
	my ($file, $hash) = @_;
	
	my $A = 1; my $C = 2; my $G = 3; my $T = 4;
	my $prev_name; my @info; my $temp;
	#slurp pwm file
	open (my $instream,  q{<}, $file) or die("Unable to open $file : $!");
		while(<$instream>){
			chomp $_;
			if(/^>/){
				$prev_name = (split/>|\s+/,$_)[1];
			}else{
				@info = split/\s+/,$_;
				if(not exists $$hash{$prev_name}){
					$$hash{$prev_name}->[0] = {(A=>$info[$A], T=>$info[$T], C=>$info[$C], G=>$info[$G])};
				}else{
					$temp = $$hash{$prev_name};
					$$hash{$prev_name}->[scalar(@$temp)] = {(A=>$info[$A], T=>$info[$T], C=>$info[$C], G=>$info[$G])};
				}
			}
		}
	close $instream;
	return;
}


sub get_pwm_id{
	my ($filename) = @_;
	my $identifier;
	my $motif_name;
	
	if($filename =~ /BVs_(.*)_(.*)_sites/){
		$motif_name = $1;
		$identifier = $2;
	}else{
		print STDERR "get_pwm_id(): ERROR - unable to recognise the input motif file name from: $filename. Aborting.\n";
		exit -1;
	}
	my $motif_id = $motif_name . '_' . $identifier;
	##heterodimers are saved by Jaspar as monomer::monomer
	##I replaced the :: with a '-' in the input file name because it's not recognised by the SGE submission	
	
	$motif_name =~ s/\-/\:\:/;
	#get the length of the motif analyzed in this iteration
	my $motif_id_postprocessed = $motif_name . '_' . $identifier;
	my $ref = $JASPAR_MOTIF{$motif_id_postprocessed};
	my $length = scalar(@$ref);	
	
	return($motif_id, $motif_id_postprocessed, $length);
}



sub write_ris_to_bed_file{
	my($in_file, $out_file, $pwm_length) = @_;
	my %hash;
	
	open (my $instream, q{<}, $in_file) or die("Unable to open $in_file : $!");
	while(<$instream>){
		chomp;
		next if($_ eq '');
		next if($_ =~ /^CHR/);
		
		my ($chr,$motif_start,$motif_end,$motif_strand,$score,$site) = (split /\t/)[0,4,5,8,9,10];
		#next if (!$chr);
		next if(  $MIN_SCORE && ($score < $MIN_SCORE) );
			
		my $pscanchip_interval_length = ($motif_end - $motif_start);
		if($pscanchip_interval_length <  $pwm_length){
			$motif_end += 1;
		}elsif($pscanchip_interval_length == $pwm_length){	
			;
		}else{
			print STDERR "ERROR: the pscanchip motif length is LARGER than the Jaspar length. Verify.\n";
			exit -1;	
		}
		my $bed_line = $chr . "\t" . $motif_start . "\t" . $motif_end . "\t" . $site . "\t" . $score;
		$hash{$bed_line} = 1;
	}
	close $instream;	
	
	open (my $outstream, q{>}, $out_file) or die("Unable to open $out_file : $!");
	foreach my $item (keys %hash){ 
		print $outstream $item, "\n"; 
	}
	close $outstream;
	return;
}


sub get_ris_intervals{
	my ($file, $pwm_length) = @_;
	my @array;
	
	open (my $instream,      q{<}, $file) or die("Unable to open $file : $!");
	while(<$instream>){
		chomp;
		next if($_ eq '');
		next if($_ =~ /^CHR/);
		
		my ($chr,$motif_start,$motif_end,$motif_strand,$score,$site) = (split /\t/)[0,4,5,8,9,10];
		#next if (!$chr);
		next if(  $MIN_SCORE && ($score < $MIN_SCORE) );
			
		my $pscanchip_interval_length = ($motif_end - $motif_start);
		if($pscanchip_interval_length <  $pwm_length){
			$motif_end += 1;
		}elsif($pscanchip_interval_length == $pwm_length){	
			;
		}else{
			print STDERR "ERROR: the pscanchip motif length is LARGER than the Jaspar length. Verify.\n";
			exit -1;	
		}
		my $bed_line = $chr . "\t" . $motif_start . "\t" . $motif_end . "\t" . $score;
		push(@array, $bed_line);
	}
	close $instream;
 	return @array;
}


# should return TRUE if the ris interval contains at least a VDR-BV; under otherwise
sub check_vdrbv_intersects_pwm_interval{
	my ($vdrbv_coords, @array) = @_;
	
	my ($vdrbv_chr,$vdrbv_pos) = split("-", $vdrbv_coords);
	#now search the array until you find an intersection
	foreach my $pwm_interval (@array){
		my ($pwm_chr, $pwm_start, $pwm_end, $score) = split("\t", $pwm_interval);
		if( ($pwm_chr eq $vdrbv_chr) && ($vdrbv_pos >= $pwm_start) && ($vdrbv_pos <= $pwm_end) ) {
			return $pwm_interval;
		}else{ 
			next;
		}	
	}
	return undef;
}

#OLD
##Get all VDR-BVs from file-----------------------------------------------------------------
#get_vdrbvs_from_file($IN_VDRBV, \%VDRBV_coords);
#my $VDRBV_initial = keys %VDRBV_coords;

#2. Fill array of pwm intervals from ris------------------------------------------------------
#my @ris_array = get_ris_intervals($RXR_VDR_PATH, $motif_length);
#3. Check the vdr-bvs against the RXR:VDR ris:------------------------------------------------
#-for each line in bed, see if vdrbv is in interval.
#If so, remove it from main hash and place in result hash
#print "Working on $full_motif_id:\n";
#foreach my $item (keys %VDRBV_coords){
#	my $best_scoring_intersecting_pwm = check_vdrbv_intersects_pwm_interval($item, @ris_array);
#	if ($best_scoring_intersecting_pwm){
#		print STDERR '.';
#		my ($chr, $start, $stop, $score) = split("\t", $best_scoring_intersecting_pwm);
#		#TODO you might want to save more than the score.
#		$RESULTS{$item}{$motif_string} = $score;
#		delete($VDRBV_coords{$item});
#	}
#}
#print "\n";
#print "Initial VDR-BVs: $VDRBV_initial\n";
#my $VDRBV_left = keys %VDRBV_coords;
#print "VDR-BV left after assigning to $full_motif_id at min score thrs: $MIN_SCORE: $VDRBV_left\n";

#old
	##2. Fill array of pwm intervals from ris------------------------------------------------------
	#my @ris_array = get_ris_intervals($RIS_FILE_PATH, $motif_length);	
#	foreach my $vdrbv (keys %VDRBV_coords){
#		my $best_scoring_intersecting_pwm = check_vdrbv_intersects_pwm_interval($vdrbv, @ris_array);
#		if ($best_scoring_intersecting_pwm){
#			print STDERR '.';
#			my ($chr, $start, $stop, $score) = split("\t", $best_scoring_intersecting_pwm);
#			$results_allpwms{$vdrbv}{$motif_string}{'SCORE'}      = $score;
#			$results_allpwms{$vdrbv}{$motif_string}{'PWM_LENGTH'} = $motif_length;
#			delete($VDRBV_coords{$vdrbv});
#		}
#	}
#	print "\n";		
#}
