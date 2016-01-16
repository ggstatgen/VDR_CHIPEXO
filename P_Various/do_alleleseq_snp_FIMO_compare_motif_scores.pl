#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#inspired by
#the mcdaniell, ewan birney paper
#http://www.sciencemag.org/content/suppl/2010/03/18/science.1184655.DC1/McDaniell.SOM.pdf

#For each allele-specific event, I want to build a graph to correlate the amount of allele specificity with the differential of the motif score log(pat_score/mat_score)
#do they correlate?
#inputs: 
#1 two text files, output of fimo, with the VDR motif instances observed around the async snps
#2 async snp file, output of alleleseq 

#The final purpose is to draw a scatterplot showing the relation of amount of "asynchronicity" in the binding VS the differential of the motif score.
#The data point is the allele-specific site. 

#ALLELESEQ FILE STRUCTURE
#chrm    snppos          ref     mat_gtyp        pat_gtyp        c_gtyp  phase   mat_all pat_all cA      cC      cG      cT      winning SymCls  SymPval BindingSite     cnv
#1       797440  T       C       T       Y       HETERO  None    None    0       0       0       8       ?       Asym    0.0078125       0       1.0
#1       1314015 C       Y       C       Y       HETERO  None    None    0       0       0       8       ?       Asym    0.0078125       0       1.0
#1       1371459 A       R       G       R       HETERO  None    None    0       0       6       0       ?       Asym    0.03125 0       1.0
#1       1509156 A       A       G       R       HETERO  None    None    0       0       5       0       ?       Sym     0.0625  0       1.0
#1       1706136 T       Y       T       Y       HETERO  None    None    0       0       0       6       ?       Asym    0.03125 0       1.0

#FIMO FILE STRUCTURE
#pattern name   sequence name   start   stop    strand  score   p-value q-value matched sequence
#0       12:26451708-26452238    76      90      +       20.844  7.85e-10        0.00047 AGGTCACTGGGTTCA
#0       19:40120439-40120620    99      113     -       20.844  7.85e-10        0.00047 AGGTCACTGGGTTCA
#0       6:143718651-143718756   70      84      +       20.844  7.85e-10        0.00047 AGGTCACTGGGTTCA
#0       19:40164228-40164391    33      47      +       20.844  7.85e-10        0.00047 AGGTCACTGGGTTCA


my $in_asymsnps;
my $in_pat_fimo;
my $in_mat_fimo;
my $SLOP;
GetOptions(
        'snps=s'     =>\$in_asymsnps,
        'pat=s'     =>\$in_pat_fimo,
        'mat=s'      =>\$in_mat_fimo
);
if(!$in_asymsnps){
	print "USAGE: do_alleleseq_snp_FIMO_compare_motif_scores.pl -snps=<ASYM_SNPS> -pat=<FIMO_PAT> -mat=<FIMO_MAT>\n";
	print "<ASYM_SNPS> file interestingHets.txt from AlleleSeq\n";
    print "<FIMO_PAT> tab separated file from FIMO listing VDR motif instances + scores, paternal\n";
    print "<FIMO_MAT> tab separated file from FIMO listing VDR motif instances + scores, maternal\n";
    exit 1;
}
if(!$in_pat_fimo){
	print "USAGE: do_alleleseq_snp_FIMO_compare_motif_scores.pl -snps=<ASYM_SNPS> -pat=<FIMO_PAT> -mat=<FIMO_MAT>\n";
	print "<ASYM_SNPS> file interestingHets.txt from AlleleSeq\n";
    print "<FIMO_PAT> tab separated file from FIMO listing VDR motif instances + scores, paternal\n";
    print "<FIMO_MAT> tab separated file from FIMO listing VDR motif instances + scores, maternal\n";
    exit 1;
}
if(!$in_mat_fimo){
	print "USAGE: do_alleleseq_snp_FIMO_compare_motif_scores.pl -snps=<ASYM_SNPS> -pat=<FIMO_PAT> -mat=<FIMO_MAT>\n";
	print "<ASYM_SNPS> file interestingHets.txt from AlleleSeq\n";
    print "<FIMO_PAT> tab separated file from FIMO listing VDR motif instances + scores, paternal\n";
    print "<FIMO_MAT> tab separated file from FIMO listing VDR motif instances + scores, maternal\n";
    exit 1;
}

#get a hash "allele-specific-site"->score for asynchronicity

#chromosome =>
#	snp_position =>
#		async/synch


open (my $instream,  q{<}, $in_asymsnps) or die("Unable to open $in_asymsnps : $!");

my %interval_to_scores; # hash of hashes which contains everything
while(<$instream>){
	chomp;
	next if($_ =~ /^chrm/);
	my ($chr,$snppos,$symCls,$symPval,$bindingSite) = (split /\t/)[0,1,14,15,16];
	next if(!$chr);
	next if(!$snppos);
	next if(!$symCls);
	next if(!$symPval);
	next if($symCls eq 'Weird');
	#next if ($symCls ne 'Asym'); #you can probably build two graphs, one for async and one for sync, showing that hopefuly the first has corr, the second doesn't
	#next if ($bindingSite ne '1');
	
	my $event =  $chr . "-" . $snppos; #the event is uniquely identified by the pair (chromosome,position)
	$interval_to_scores{$chr}{$snppos}{$symCls} = $symPval;
}
close $instream;

#get a hash "allele specific event"->fimo score


		


my $counter = 0;
open ($instream,  q{<}, $in_pat_fimo) or die("Unable to open $in_pat_fimo : $!");
while(<$instream>){
	next if($_ =~ /^\#/);
	chomp;
	my ($coord,$score) = (split /\t/)[1,5];
	next if(!$coord);
	next if(!$score);
	
	my ($chr,$interval) = split(':',$coord);
	my ($start,$stop) = split('-', $interval);
	
	foreach my $this_chr (keys %interval_to_scores){
		if($chr eq $this_chr){
			for my $this_snp (keys %{ $interval_to_scores{$this_chr} }){
				if( ($this_snp >= $start) && ($this_snp <= $stop) ){
					$interval_to_scores{$this_chr}{$this_snp}{"pat"} = $score;
					$counter++;
				}
			}
		}
	}
}
close $instream;

#same with maternal
open ($instream,  q{<}, $in_mat_fimo) or die("Unable to open $in_mat_fimo : $!");
while(<$instream>){
	next if($_ =~ /^\#/);
	chomp;
	my ($coord,$score) = (split /\t/)[1,5];
	next if(!$coord);
	next if(!$score);
	
	my ($chr,$interval) = split(':',$coord);
	my ($start,$stop) = split('-', $interval);
	
	foreach my $this_chr (keys %interval_to_scores){
		if($chr eq $this_chr){
			for my $this_snp (keys %{ $interval_to_scores{$this_chr} }){
				if( ($this_snp >= $start) && ($this_snp <= $stop) ){
					$interval_to_scores{$this_chr}{$this_snp}{"mat"} = $score;
				}
			}
		}
	}
}
close $instream;

#create outfile for matlab?
#chromosome,	position, event, event_pval,	mat_score, pat_score    
print "chrm\tposition\tevent\tevent_pval\tpat_score\tmat_score\tlogdiff_score\n";
foreach my $this_chr ( sort keys %interval_to_scores ) {
    foreach my $this_snp ( sort keys %{ $interval_to_scores{$this_chr} } ) {
		foreach my $event ( sort keys  %{ $interval_to_scores{$this_chr}{$this_snp} }){
			print $this_chr . "\t" . 
			      $this_snp . "\t" .
			      $event    . "\t" . 
			      -log($interval_to_scores{$this_chr}{$this_snp}{$event}) . "\t" . 
			      #$interval_to_scores{$this_chr}{$this_snp}{"pat"}  . "\t" .
			      #$interval_to_scores{$this_chr}{$this_snp}{"mat"}  . "\t" . 
			      sprintf('%.5f', log($interval_to_scores{$this_chr}{$this_snp}{"pat"}/$interval_to_scores{$this_chr}{$this_snp}{"mat"})) . "\n"
			      if(($event =~ /sym/i) && $interval_to_scores{$this_chr}{$this_snp}{"mat"});
		}
    }
}
