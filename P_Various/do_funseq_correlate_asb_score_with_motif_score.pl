#!/usr/bin/perl
use strict;
use warnings;
#use File::Basename;
#use Getopt::Long;
use List::Util qw(min max);

#8/1/2015
#CPP comment
#I was expecting an analysis comparing the extent of binding change with the extent of matching to the pwm VDR:RXR motif

#THIS SCRIPT IS UNFINISHED - there will not be a correlation for most points, those falling into less interesting nucleotides. Discuss with Chris first.


#So here I want to collect data pairs of the kind (genotype_change, phenotype_change) and then make a scatter plot
#genotype_change is the change in PWM score for motif break events
#phenotype_change is the change in read depth for that event

#also check this image
#the mcdaniell, ewan birney paper
#http://www.sciencemag.org/content/suppl/2010/03/18/science.1184655.DC1/McDaniell.SOM.pdf

#For each allele-specific event, I want to build a graph to correlate the amount of allele specificity (p-value?) with the differential of the motif score log(pat_score/mat_score)

#there will be events called as allele-specific at the same position for multiple samples. These will have different p-values for each sample. Which one to use? Use the best or the worst? ---> currently use the smalles


#do they correlate?
#inputs: 
#1 - collective output of alleleseq, using all samples. Doesn't matter that samples are not labelled. You will build a hash with "coord" -> {p-values} and use the smallest p-value
#2 - output.vcf with motif break scores out of funseq


my $infile_ASB_collective = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/ENCODE_FILTERED/ASB_EVENTS_BY_SAMPLE/interestingHets_collective_EBLfiltered.txt";
#format
#1       1080920 G       S       W       R       PHASED  G       A       2       0       7       0       M       Sym     0.1796875       1       1.0
#1       1080925 G       C       S       S       PHASED  C       G       0       0       8       0       P       Asym    0.0078125       1       1.0

my $infile_funseq = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/Output.vcf";

my $motif_name = "VDR_JASPAR";

#read in p-values for asym events and save in hash 
my %asbevent_to_pvalues_asym;
my %asbevent_to_pvalues_sym;
open (my $instream,  q{<}, $infile_ASB_collective) or die("Unable to open $infile_ASB_collective : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^chrm/);
	
	my ($chr,$snppos,$symCls,$symPval) = (split /\t/)[0,1,14,15];
	next if(!$chr);
	next if(!$snppos);
	next if(!$symCls);
	next if(!$symPval);
	next if($symCls eq 'Weird');
	
	next if ($symCls ne 'Asym'); #you can probably build two graphs, one for async and one for sync, showing that hopefuly the first has corr, the second doesn't
	#next if ($bindingSite ne '1');
	
	my $event =  $chr . '-' . $snppos; #the event is uniquely identified by the pair (chromosome,position)
	
	if($symCls eq 'Asym'){
		if(!$asbevent_to_pvalues_asym{$event}){
			$asbevent_to_pvalues_asym{$event} = $symPval;
		}else{
			$asbevent_to_pvalues_asym{$event} = $asbevent_to_pvalues_asym{$event} . ',' . $symPval;
		}		
	}elsif($symCls eq 'Sym'){
		if(!$asbevent_to_pvalues_sym{$event}){
			$asbevent_to_pvalues_sym{$event} = $symPval;
		}else{
			$asbevent_to_pvalues_sym{$event} = $asbevent_to_pvalues_sym{$event} . ',' . $symPval;
		}				
	}

}
close $instream;



#grab each motif breaking event for VDR_JASPAR - build the event; output pair p-value / log(paternal score/maternal score)
open ($instream,  q{<}, $infile_funseq) or die("Unable to open $infile_funseq : $!");
while(<$instream>){
	chomp;
	my $info_motifbr;
	
	next if($_ =~ /^\#/);
	next if($_ eq '');
	next unless($_ =~ /MOTIF/);
	
	#get info field
	my ($chr, $pos, $info) = (split /\t/)[0,1,7];
	$chr =~ s/chr(.*)/$1/;
	
	my $event = $chr . '-' . $pos;
	
	my @info = split(";", $info);
	foreach my $item (@info){
		$info_motifbr = $item if($item =~ /^MOTIFBR/);
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
			#next unless( $TF_motif eq $motif_name);
		
			#here you know you have a break in the VDR motif/ Build the log score
			#my $break_score = log ( ($TF_ref_score + 1) / ($TF_alt_score + 1) ); 
			
			my $break_score = $TF_alt_score - $TF_ref_score; 
			
			#associate this to event
			my $pvalues = $asbevent_to_pvalues_asym{$event};
			my $candidate_pval = max split(",", $pvalues);
			print $candidate_pval . "\t" . $break_score . "\n";
		}
	}
}
close $instream;


sub log2 {
	my $n = shift;
	return log($n)/log(2);
}
