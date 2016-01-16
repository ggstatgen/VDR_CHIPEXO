#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#I want to create a bar plot for the last figure of the paper where I list one bar per category
#dhs
#enhancer
#hot region
#motif disrupt any
#motif disrupt vdr
#recur
#gwas
#PICS
#etc

#19/jan/2015 chris p is confused by the value in VDR peak break. Should I only list VDR_JASPAR breaks? One column for DREME and one column for JASPAR?

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
#$infile = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/Output_noDBRECUR.vcf"; #I have removed all the DBRECUR instances not to get confused
#$FILTER_RECUR = 1;
#my $infile_pics = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/BROAD_PICS/BROAD_candidate_causal_snps_39immune_nonimmune_diseases_plus_enh_annot_masterfile9_hg19.vcf.gz";
#my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools2-2.20.1/bin/bedtools";

if(!$infile){
	print "USAGE: do_funseq_produce_barplot.pl -i=<INFILE> -vdr -recur\n";
    print "<INFILE> vcf output of funseq2\n";
    print "(optional)<vdr> whether to only use ASB variants in a VDR peak (default=no)\n";
    print "(optional)<recur> whether to only use ASB variants which recur (default=no)\n";
    exit 1;
}


my $counter_all = 0; # this is either a) all asb annotated variants b) all asb annotated variants in VDR peaks c) all asb annotated variants which recur or b)+c)
#NCENC
my $counter_ncenc_enhancer = 0;
my $counter_ncenc_dhs = 0;
my $counter_ncenc_lncrna = 0;
my $counter_ncenc_tfp = 0;
my $counter_ncenc_tfp_vdr = 0;
my $counter_ncenc_tfm = 0; 
my $counter_ncenc_tfm_vdr = 0;
my $counter_ncenc_pseudogene = 0;
#BOOLS
my $counter_hot = 0;
my $counter_recur = 0;
my $counter_ppi_hub = 0;
my $counter_sen = 0;
my $counter_usen = 0;
my $counter_ucon = 0;
#MOTIFBR
my $counter_motif_break = 0;
my $counter_motif_break_vdr = 0;
#GENE
my $counter_gene_intron = 0; #GENE=CASP9(Intron)
my $counter_gene_distal = 0; #GENE=MFAP2(Distal)
my $counter_gene_promoter = 0; #GENE=MFAP2(Promoter) GENE=UBR4(Intron&Promoter)
my $counter_gene_utr = 0; #GENE=MFAP2(UTR)
my $counter_gene_medial = 0;
#not annotated by any of the above
my $counter_none = 0;

my %id_to_info;
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/);
	next if($_ eq '');
	
	#get info field
	my ($chr, $pos, $rsID, $info) = (split /\t/)[0,1,2,7];
	#filters
	if($FILTER_RECUR){
		next unless($info =~ /RECUR/);
	}
	if($FILTER_VDR_ONLY){
		next unless($info =~ /NCENC/);
	}
	my $id = $chr . '-' . $pos; #recurrent is  printed once per sample in output.vcf
	if($id_to_info{$id}){
		next;
	}
	$id_to_info{$id} = $info;
}

foreach my $id (keys %id_to_info){
	my $info = $id_to_info{$id};
	my $info_motifbr;
	my $info_motifg;
	my $info_ncenc;
	my $info_gene;
	my $info_recur;
	
	my $no_gene; my $no_ncenc; my $no_motifbr; my $no_bool;
	
	#non booleans: ncenc, gene, motifbr, motifg
	my @info = split(";", $info);
	foreach my $item (@info){
		$info_ncenc = $item if($item =~ /^NCENC/);
		$info_gene = $item if($item =~ /^GENE/);
		$info_motifbr = $item if($item =~ /^MOTIFBR/);
		$info_motifg = $item if($item =~ /^MOTIFG/);
		$info_recur = $item if($item =~ /^RECUR/);
	}
	#VDR peak only filter 
	if($FILTER_VDR_ONLY){
		next unless($info_ncenc =~ /TFP(.*)VDR/);
	}
	$counter_all++;

	#see here http://info.gersteinlab.org/Funseq2#F._Output_files
	#for formats
	#if ncenc exists, and if it contains TFP and TFM, I want to separate the VDR from the rest. In other words, I don't want 
	#to count the VDR twice
	
	if($info_ncenc){	
		$counter_ncenc_dhs++ if( $info_ncenc =~ /DHS/);
		$counter_ncenc_enhancer++ if( $info_ncenc =~ /Enhancer/);
		$counter_ncenc_lncrna++ if( $info_ncenc =~ /ncRNA/);
		$counter_ncenc_pseudogene++ if( $info_ncenc =~ /Pseudogene/);
		
		if( $info_ncenc =~ /TFP/){
			my $SEEN_VDR = 0; my $SEEN_ENCODE_TF = 0;
			
			my($tag,$info) = split('=',$info_ncenc);
			my @fields = split(',',$info);
			foreach my $item (@fields){
				if($item =~ /TFP\((.+)\|/){
					if($1 eq 'VDR'){
						$SEEN_VDR = 1;
					}else{
						$SEEN_ENCODE_TF = 1;
					}
				}
			}	
			$counter_ncenc_tfp++ if($SEEN_ENCODE_TF);
			$counter_ncenc_tfp_vdr++ if($SEEN_VDR);
		}
		
		if( $info_ncenc =~ /TFM/){
			my $SEEN_VDR_JASPAR = 0; my $SEEN_ENCODE_TF = 0;
			
			my($tag,$info) = split('=',$info_ncenc);
			my @fields = split(',',$info);
			foreach my $item (@fields){
				if($item =~ /^TFM\((\w+)\|/){
					#print $1 . "\n";
					if($1 eq 'VDR'){
						#do subfields here if you want to count other vdr custom motifs
						#$SEEN_VDR_JASPAR = 1 if($item =~ /TFM\(VDR\|VDR_JASPAR/);
						$SEEN_VDR_JASPAR = 1 if($item =~ /JASPAR/);
					}else{
						$SEEN_ENCODE_TF = 1;
					}
				}
			}	
			$counter_ncenc_tfm++ if($SEEN_ENCODE_TF);
			$counter_ncenc_tfm_vdr++ if($SEEN_VDR_JASPAR);
		}
	}else{
		$no_ncenc = 1;
	}
	
	#19/1/2015
	#Should I only count instances of VDR_JASPAR break for coherence with FIGURE 3?
	
	if($info_motifbr){
		#MOTIFBR=THAP1#THAP1_disc1_8mer#7680676#7680701#+#21#0.020000#0.680000,VDR#VDR_XXmotif#7680691#7680706#+#6#0.03084#0.89236
		my $SEEN_VDR_JASPAR = 0; my $SEEN_ENCODE_TF = 0;
		
		my($tag,$info) = split('=',$info_motifbr);
		my @fields = split(',',$info);	
		foreach my $item (@fields){
			if($item =~ /^VDR/){
				$SEEN_VDR_JASPAR = 1 if($item =~ /JASPAR/);
			}else{
				$SEEN_ENCODE_TF = 1;
			}
		}	
		$counter_motif_break++ if($SEEN_ENCODE_TF);
		$counter_motif_break_vdr++ if($SEEN_VDR_JASPAR);	
	}else{
		$no_motifbr = 1;
	}
	
	
	#For noncoding variants, intron, promoter, UTR, Distal and Medial tags are annotated.
	if($info_gene){
		$counter_gene_promoter++ if($info_gene =~ /GENE(.*)Promoter/);
		$counter_gene_intron++ if($info_gene =~ /GENE(.*)Intron/);
		$counter_gene_distal++ if($info_gene =~ /GENE(.*)Distal/);
		$counter_gene_medial++ if($info_gene =~ /GENE(.*)Medial/);
		$counter_gene_utr++ if($info_gene =~ /GENE(.*)UTR/);
	}else{
		$no_gene = 1;
	}
	

	
	#booleans
	#sensitive, ultrasensitive, ultraconserved, recur, hot, hub
	my $at_least_one_bool;
	if($info =~ /HOT/){
		$counter_hot++;
		$at_least_one_bool = 1;
	}
	if($info =~ /USEN/){
		$counter_usen++;
		$at_least_one_bool = 1;
	}
	if($info =~ /UCON/){
		$counter_ucon++;
		$at_least_one_bool = 1;
	}
	if( ($info =~ /^SEN/) ||  ($info =~ /;SEN/) ){
		$counter_sen++;
		$at_least_one_bool = 1;
	}
	if($info =~ /RECUR/){
		$counter_recur++;
		$at_least_one_bool = 1;
	}
	if($info =~ /HUB/){
		$counter_ppi_hub++;
		$at_least_one_bool = 1;
	}
	$no_bool = 1 unless($at_least_one_bool);
	
	if($FILTER_VDR_ONLY){
		$counter_none++ if($no_bool && $no_gene && $no_motifbr);
	}elsif($FILTER_VDR_ONLY && $FILTER_RECUR){
		$counter_none++ if($no_gene && $no_motifbr);
	}elsif($FILTER_RECUR){
		$counter_none++ if($no_gene && $no_motifbr && $no_ncenc);
	}else{
		$counter_none++ if($no_bool && $no_gene && $no_motifbr && $no_ncenc);
	}
}

#fractions
#$counter_ncenc_dhs = $counter_ncenc_dhs / $counter_all ;
#$counter_ncenc_enhancer = $counter_ncenc_enhancer / $counter_all;
#$counter_ncenc_lncrna = $counter_ncenc_lncrna / $counter_all;
#$counter_ncenc_tfp = $counter_ncenc_tfp / $counter_all;
#$counter_ncenc_tfm = $counter_ncenc_tfm / $counter_all;
#$counter_ncenc_tfp_vdr = $counter_ncenc_tfp_vdr / $counter_all;
#$counter_ncenc_tfm_vdr = $counter_ncenc_tfm_vdr / $counter_all;
#
#$counter_gene_promoter = $counter_gene_promoter / $counter_all;
#$counter_gene_intron = $counter_gene_intron / $counter_all;
#$counter_gene_distal = $counter_gene_distal / $counter_all;
#$counter_gene_medial = $counter_gene_medial / $counter_all;
#$counter_gene_utr = $counter_gene_utr / $counter_all;
#
#$counter_motif_break = $counter_motif_break / $counter_all;
#$counter_motif_break_vdr = $counter_motif_break_vdr / $counter_all;
#
#$counter_hot = $counter_hot / $counter_all;
#$counter_sen = $counter_sen / $counter_all;
#$counter_usen = $counter_usen / $counter_all;
#$counter_ucon = $counter_ucon / $counter_all;
#$counter_ppi_hub = $counter_ppi_hub / $counter_all;
#$counter_recur = $counter_recur / $counter_all;
#$counter_none = $counter_none / $counter_all;

#header
print "ANNOTATION\tOBSERVED\n";
print 'ALL_ANNOTATED_VDR_VARIANTS' . "\t" . $counter_all . "\n";
print 'In_DHS' . "\t" . $counter_ncenc_dhs . "\n";
print 'In_Enhancer' . "\t" . $counter_ncenc_enhancer . "\n";
print 'In_lncRNA' . "\t" . $counter_ncenc_lncrna . "\n";
print 'In_ENCODE_TF_Peak' . "\t" . $counter_ncenc_tfp . "\n";
print 'In_ENCODE_TF_Motif' . "\t" . $counter_ncenc_tfm . "\n";
print 'ENCODE_TF_Motif_Break' . "\t" . $counter_motif_break . "\n";
print 'In_VDR_Peak' . "\t" . $counter_ncenc_tfp_vdr . "\n";
print 'In_VDR:RXR_Motif' . "\t" . $counter_ncenc_tfm_vdr . "\n";
print 'VDR:RXR_Motif_Break' . "\t" . $counter_motif_break_vdr . "\n";
print 'In_Gene_Promoter' . "\t" . $counter_gene_promoter . "\n";
print 'In_Gene_Intron' . "\t" . $counter_gene_intron . "\n";
print 'In_Gene_Distal' . "\t" . $counter_gene_distal . "\n";
print 'In_Gene_Medial' . "\t" . $counter_gene_medial . "\n";
print 'In_Gene_UTR' . "\t" . $counter_gene_utr . "\n";
print 'In_Hot_region' . "\t" . $counter_hot . "\n";
print 'Is_Sensitive' . "\t" . $counter_sen . "\n";
print 'Is_Ultrasensitive' . "\t" . $counter_usen . "\n";
print 'Is_Ultraconserved' . "\t" . $counter_ucon . "\n";
print 'Is_Recurrent' . "\t" . $counter_recur . "\n";
print 'None_of_the_above' . "\t" . $counter_none . "\n";