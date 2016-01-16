#!/usr/bin/perl
use strict;
use warnings;


#get a latex table summary from the funseq/ broad intersection results

#bedtools intersect -wo -a Output_recur_allfields_s.vcf.gz -b /net/isi-scratch/giuseppe/VDR/ALLELESEQ/BROAD_PICS/BROAD_candidate_causal_snps_39immune_nonimmune_diseases_plus_enh_annot_masterfile9_hg19.vcf.gz > Output_recur_INTERSECT_BROAD_PICS.tsv


#my $INPUT = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/Output_sorted_INTERSECT_BROAD_PICS.tsv";
#my $INPUT = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/Output_recur_INTERSECT_BROAD_PICS.tsv";
my $INPUT = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/vdrbv_pics/Output_sorted_INTERSECT_BROAD_PICS.tsv";

#format
#chr1    2539555 rs4278312       A       G       100     PASS    SAMPLE=interestingHets_NA19213_EBLfiltered_hg19;GERP=-0.859;CDS=No;HUB=MMEL1:PPI(0.142);GENE=MMEL1(Intron);NCDS=0.0625421843697522      chr1    2539555 rs4278312       - -       0.0505  Celiac_disease  RISK_ALLELE=T;ANNOTATION=none;NEAREST_GENE=MMEL1;EQTL=none;EQTL_DIR=NA;TOP_ENH=none     1

#get fields 0,1,2,3,4,7,13,14,15
#further process fields 7 and 15 into columns

my $header = "CHR\tPOS\tRS_ID\tFUNSEQ_ANN\tFUNSEQ_SCORE\tPICS_SCORE\tPICS_TRAIT\tPICS_ANN";

print $header, "\n";
open (my $instream,  q{<}, $INPUT) or die("Unable to open $INPUT : $!");
while(<$instream>){
	chomp;
	my ($chr,$snppos,$rsid, $ref, $alt, $funseq_info, $pics_score, $pics_trait, $pics_info) = (split /\t/)[0,1,2,3,4,7,13,14,15];
	
	#further process $funseq_info and $pics_info
	my @funseq_info = split(";", $funseq_info);
	my @funseq_info_processed;
	my $funseq_score;
	foreach my $item (@funseq_info){
		if($item =~ /^SAMPLE/){ next; }
		if($item =~ /^CDS/){ next; }
		elsif($item =~ /NCDS/){
			my ($name, $temp) = split("=", $item); 
			my $value = sprintf "%.3f", $temp;
			$funseq_score = $value;
			#$funseq_score = $item;
		}
		elsif($item =~ /GERP/){ push (@funseq_info_processed, $item); }
		elsif($item =~ /NCENC/) { #quite long, simplify
			my ($name, $temp) = split("=", $item);
			my @temp = split(",", $temp);
			my %anno;my @anno;
			foreach my $anno (@temp){
				if($anno =~ /(.*)\(.*\)/){
					$anno{$1} = 1;
				}else{
					print "Warning: could not regex annotation: $anno\n";
				}
			}
			foreach my $anno (sort keys %anno){ push(@anno, $anno); }
			my $anno = join(",", @anno);
			my $new_string = $name . "=" . $anno;
			push (@funseq_info_processed, $new_string); 
		}elsif($item =~ /MOTIFBR/) { push (@funseq_info_processed, 'MOTIFBR=Y'); }
		elsif($item =~ /MOTIFG/) { push (@funseq_info_processed, 'MOTIFG=Y'); }
		elsif($item =~ /SEN/) { push (@funseq_info_processed, $item); }
		elsif($item =~ /CONS/) { push (@funseq_info_processed, $item); }
		elsif($item =~ /GENE/) { #quite long, simplify
			my ($name, $temp) = split("=", $item);
			my @temp = split(",", $temp);
			my %anno;my @anno;
			foreach my $anno (@temp){
				if($anno =~ /(.*)\(.*\)\[.*\]/){
					$anno{$1} = 1;
				}elsif($anno =~ /(.*)\(.*\)/){
					$anno{$1} = 1;
				}
				else{
					print "Warning: could not regex annotation: $anno\n";
				}
			}
			foreach my $anno (sort keys %anno){ push(@anno, $anno); }
			my $anno = join(",", @anno);
			my $new_string = $name . "=" . $anno;
			push (@funseq_info_processed, $new_string); 
		}elsif($item =~ /HUB/) { #quite long, simplify
			my @temp = split(":",$item);
			push (@funseq_info_processed, $temp[0]); 
		}elsif($item =~ /HOT/) { push (@funseq_info_processed, 'HOT=Y'); }
		elsif($item =~ /RECUR/) { push (@funseq_info_processed, 'RECUR=Y'); }
		else{ push (@funseq_info_processed, $item); }
	}
	
	my $funseq_info_processed = join(';', @funseq_info_processed);
	print $chr . "\t" . 
			$snppos . "\t" . 
			$rsid   . "\t" .
			$funseq_info_processed . "\t" . 
			$funseq_score . "\t" . 
			$pics_score . "\t" . 
			$pics_trait . "\t" . 
			$pics_info . "\n"; 
}
