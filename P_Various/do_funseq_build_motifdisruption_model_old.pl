#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;


#I want to build a histogram counting and binning how the motif disrupting events fall across the VDR:RXR motif.
#The output would be an image showing which bases seem to be disrupted more often

#Input:
#Output.vcf file from funseq
#you need two subsections of the info field, MOTIFBR and MOTIFG

#eg
#MOTIFBR=MAX#Myc_known9_8mer#102248644#102248656#-#9#0.068966#0.931034
#TF name # motif name # motif start # motif end # motif strand # mutation position # alternative allele frequency in PWM # reference allele frequency in PWM
#MOTIFG=GATA_known5#75658824#75658829#-#1#4.839#4.181
#motif name # motif start # motif end # motif strand # mutation position # sequence score with alternative allele # sequence score with reference allele. (0-based, end exclusive)


#I also filter based on positive/negative GERP scores. The definition is here  http://ucscbrowser.genap.ca/cgi-bin/hgTables?db=hg19&hgta_group=compGeno&hgta_track=allHg19RS_BW&hgta_table=allHg19RS_BW&hgta_doSchema=describe+table+schema

#Sites are scored independently. Positive scores represent a substitution deficit (i.e., fewer substitutions than the average neutral site) and thus indicate that a site may be under evolutionary constraint. Negative scores indicate that a site is probably evolving neutrally; negative scores should not be interpreted as evidence of accelerated rates of evolution because of too many strong confounders, such as alignment uncertainty or rate variance. Positive scores scale with the level of constraint, such that the greater the score, the greater the level of evolutionary constraint inferred to be acting on that site.

#We applied GERP, as implemented in the GERP++ software package, to quantify the level of evolutionary constraint acting on each site in hg19, based on an alignment of 35 mammals to hg19 with a maximum phylogenetic scope of 6.18 substitutions per neutral site. Gaps in the alignment are treated as missing data, which means that the number of substitutions per neutral site will be less than 6.18 in sites where one or more species has a gap. Thus, RS scores range from a maximum of 6.18 down to a below-zero minimum, which we cap at -12.36. RS scores will vary with alignment depth and level of sequence conservation. A score of 0 indicates that the alignment was too shallow at that position to get a meaningful estimate of constraint. Should classification into "constrained" and "unconstrained" sites be desired, a threshold may be chosen above which sites are considered "constrained". In practice, we find that a RS score threshold of 2 provides high sensitivity while still strongly enriching for truly constrained sites. 


my $infile;
my $motif_name;
my $gerp_filter;
my $ENH_ONLY;
GetOptions(
        'i=s'        =>\$infile,
        'motif=s'    =>\$motif_name,
        'gerp=s'     =>\$gerp_filter,
        'enh'   	 =>\$ENH_ONLY
);

#$infile = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/Output.vcf";
#$motif_name = "VDR_dreme";
#$ENH_ONLY = 1;
#$gerp_filter = 'pos';

if(!$infile){
	print "USAGE: do_funseq_build_motifdisruption_model.pl -i=<INFILE> -motif=<MOTIF_NAME> -gerp=<pos|neg> -enh\n";
    print "<INFILE> vcf output of funseq2\n";
    print "<MOTIF_NAME> motif name to scan for. One of [VDR_JASPAR|VDR_XXmotif|VDR_dreme]\n";
    print "(optional)<gerp> whether to restric the snp labelled with positive or negative GERP score. ONE of [pos|neg]\n";
    print "(optional)<enh> whether to only use snps with NCENC=enhancer (default=no)\n";
    print "note: gathering counts for [VDR_JASPAR|VDR_XXmotif]";
    exit 1;
}
if(!$motif_name){
	print "USAGE: do_funseq_build_motifdisruption_model.pl -i=<INFILE> -motif=<MOTIF_NAME> -gerp=<pos|neg> -enh\n";
    print "<INFILE> vcf output of funseq2\n";
    print "<MOTIF_NAME> motif name to scan for. One of [VDR_JASPAR|VDR_XXmotif|VDR_dreme]\n";
    print "(optional)<gerp> whether to restric the snp labelled with positive or negative GERP score. ONE of [pos|neg]\n";
    print "(optional)<enh> whether to only use snps with NCENC=enhancer (default=no)\n";
    print "note: gathering counts for [VDR_JASPAR|VDR_XXmotif]";
    exit 1;
}

my %motifbr_pos_to_frequency;
my %motifg_pos_to_frequency;
#my @motif_break_positions;
#my @motif_gain_positions;
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	my $info_motifbr;
	my $info_motifg;
	my $info_ncenc;
	my $info_gerp;
	
	next if($_ =~ /^\#/);
	next if($_ eq '');
	next unless($_ =~ /MOTIF/);
	
	if($ENH_ONLY){ next unless($_ =~ /NCENC/); }
	
	#get info field
	my ($chr, $pos, $info) = (split /\t/)[0,1,7];
	my $ID = $chr . ' ' . $pos;
	
	my @info = split(";", $info);
	foreach my $item (@info){
		$info_motifbr = $item if($item =~ /^MOTIFBR/);
		$info_motifg = $item if($item =~ /^MOTIFG/);
		$info_ncenc = $item if($item =~ /^NCENC/);
		$info_gerp = $item if($item =~ /^GERP/);
	}
	
	#EHNANCER FILTER
	if($ENH_ONLY){ next unless($info_ncenc =~ /Enhancer/); }
	#GERP FILTER
	if($gerp_filter){
		my ($gerp_tag,$gerp_score) = split("=", $info_gerp);		
		#some snps have a gerp score set to '.'. I will assign zero to them so they will be skipped in the control down here
		if($gerp_score eq '.'){ $gerp_score = 0; }
		
		if($gerp_filter eq 'pos'){
			next unless($gerp_score > 0);
		}elsif($gerp_filter eq 'neg'){
			next unless($gerp_score < 0);	
		}else{
			print "Error: gerp filter option: $gerp_filter not recognised. Aborting. \n";
			exit -1;
		}
	}
	
	
	#MOTIFBR=VDR#VDR_JASPAR#139340766#139340781#+#13#0.000#0.900,VDR#VDR_XXmotif#139340766#139340781#+#13#0.04338#0.80145
	if($info_motifbr){
		my ($motif_change_type,$motif_change_byTF) = split("=", $info_motifbr);
		my @motif_change_byTF = split(",", $motif_change_byTF);
		foreach my $item (@motif_change_byTF){
			my ($TF,
				$TF_motif,
				$TF_motif_start,
				$TF_motif_end,
				$TF_motif_strand,
				$TF_snp_position,
				$TF_alt_score,
				$TF_ref_score) =  split("#", $item);
			next unless( $TF_motif eq $motif_name);
			my $motif_size = $TF_motif_end - $TF_motif_start;
			
			if($TF_motif_strand =~ /-/){
				my $TF_snp_position_RC = ($motif_size - $TF_snp_position + 1 );
				#push(@motif_break_positions,  $TF_snp_position_RC );
				my $break_position = 'pos_' . $TF_snp_position_RC;
				$motifbr_pos_to_frequency{$break_position} += 1;
			}else{
				#push(@motif_break_positions, $TF_snp_position);
				my $break_position = 'pos_' . $TF_snp_position;
				$motifbr_pos_to_frequency{$break_position} += 1;
			}
		}
	}
	
	#MOTIFG=GATA_known5#75658824#75658829#-#1#4.839#4.181
	if($info_motifg){
		my ($motif_change_type,$motif_change_byTF) = split("=", $info_motifg);
		my @motif_change_byTF = split(",", $motif_change_byTF);
		foreach my $item (@motif_change_byTF){
			my ($TF_motif,
				$TF_motif_start,
				$TF_motif_end,
				$TF_motif_strand,
				$TF_snp_position,
				$TF_alt_score,
				$TF_ref_score) =  split("#", $item);
			next unless( $TF_motif eq $motif_name);
			my $motif_size = $TF_motif_end - $TF_motif_start;
			
			if($TF_motif_strand =~ /-/){
				my $TF_snp_position_RC = ($motif_size - $TF_snp_position + 1 );
				#push(@motif_gain_positions,  $TF_snp_position_RC );
				my $gain_position = 'pos_' . $TF_snp_position_RC;
				$motifg_pos_to_frequency{$gain_position} += 1;
			}else{
				#push(@motif_gain_positions, $TF_snp_position);
				my $gain_position = 'pos_' . $TF_snp_position;
				$motifg_pos_to_frequency{$gain_position} += 1;
			}		
		}
	}
}
close $instream;

print "MOTIFBR_POSITION\tMOTIFBR_FREQ\n";
foreach my $item (sort  keys %motifbr_pos_to_frequency){
	print $item . "\t" . $motifbr_pos_to_frequency{$item} . "\n";
}
print "MOTIFG_POSITION\tMOTIFG_FREQ\n";
foreach my $item (sort  keys %motifg_pos_to_frequency){
	print $item . "\t" . $motifg_pos_to_frequency{$item} . "\n";	
}


#print "VDR_JASPAR_MOTIFBR\n";
#foreach my $item (@motif_break_positions){
#	print $item, "\n";
#}
#
#print "VDR_JASPAR_MOTIFG\n";
#foreach my $item (@motif_gain_positions){
#	print $item, "\n";
#}

