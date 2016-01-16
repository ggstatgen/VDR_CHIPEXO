#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#can I melt this with the other?

#21/1/2015 
#I want to collect all GERP scores for all variants in VDR_JASPAR motifs (or subset of theses)
#Then, I want to plot a boxplot for each of the motif positions with the GERP scores
#label GERP scores for variants which break the motif


#Output is as follows:
#type	position	gerp_score
#break
#no_break (lob? gob? you'd need the phenotype for this)

#Input:
#1)######################################################
#Output.vcf file from funseq - I need lines tagged with GERP(all?), NCENC and TFM VDR_JASPAR (labelled as falling under the JASPAR motif)

#2)####################################################
#The motif file I originally used for the funseq analysis. I need these to see if the motif is on the + or - strand, in order to compute WHERE
#in this motif the SNP is hitting (position)

my $infile;
my $INPUT_MOTIF_FILE = '/net/isi-scratch/giuseppe/tools/funseq2-1.0/data/ENCODE.tf.bound.union.bed';
my $motif_name = 'VDR_JASPAR';
my $ENH_ONLY;
my $subset_bed; #if you only want to collect data for positions in an external bed file (eg CLASS I sites)
my $sym_variants; #flag, set it if you're dealing with sym variants

#vector of possible motif positions - this will be hardcoded to a 15nt motif
#I need this to insert 'NA's in the R output when for a position there is no data at all
my @MOTIF_POSITIONS = ("pos_1", "pos_2", "pos_3", "pos_4", "pos_5", "pos_6", "pos_7", "pos_8", "pos_9", "pos_10", "pos_11", "pos_12", "pos_13", "pos_14", "pos_15");
my $MOTIF_SIZE = 15;
GetOptions(
        'i=s'      =>\$infile,
        'subset=s' =>\$subset_bed,
        'sym'      =>\$sym_variants,
        'enh'      =>\$ENH_ONLY
); 
$infile = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/Output.vcf";
#$infile = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_sym_allsamples_ancestral/Output.vcf";
#$ENH_ONLY = 1;
$subset_bed = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_motif_scan_pscanchip/Pscanchip_occurrences_VDRRXR_jaspar_classI.bed";
#$sym_variants = 1;

if(!$infile){
	print "USAGE: do_funseq_collect_GERP.pl -i=<INFILE> -subset=<SUBSET_BED> -sym -enh\n";
    print "<INFILE> vcf output of funseq2\n";
    print "(optional)<SUBSET_BED> only collect data for VDR-BVs breaking motifs intersecting with this bed file(default=no)\n";
    print "(optional)<sym> flag; set it if you are dealing with SYM variants from alleleseq (default=no)";
    print "(optional)<enh> whether to only use snps with NCENC=enhancer (default=no)\n";
    exit 1;
}
my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
my $output;

if($subset_bed){
	$output  = $directory . $basename . '_GERPscores_subsetbed.Rdata';
}elsif($ENH_ONLY){
	$output  = $directory . $basename . '_GERPscores_enhonly.Rdata';
}else{
	$output  = $directory . $basename . '_GERPscores_all.Rdata';
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
#1 get all the motif instance coordinates with their orientation.
#####################
#This is needed to find the position where the SNP hits in the motif sequence
my %vdr_motif_to_strand;
open (my $instream,  q{<}, $INPUT_MOTIF_FILE) or die("Unable to open $INPUT_MOTIF_FILE : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	
	my @fields = split("\t", $_);
	#format is as follows
	#chr1    1310676 1310683 VDR_dreme       .       +       VDR
	#chr1    1310714 1310729 VDR_XXmotif     .       +       VDR
	#chr1    1310747 1310762 VDR_JASPAR      .       +       VDR
	next unless($fields[3] eq $motif_name);
	my $identifier = $fields[0] . ':' . $fields[1] . '-' . $fields[2];     #chr7:5013558-5013573
	$vdr_motif_to_strand{$identifier} = $fields[5];
}
close $instream;

#########################
#2 pick alleleseq snps hitting motifs from funseq
#get their GERP score
#########################
my %motifmodel_motifpos2geomicpos2GERP;
my %motifmodel_motifpos2geomicpos2GERP_MOTIFBR;
open ($instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/);
	next if($_ eq '');
	my $info_ncenc;
	my $info_motifbr;
	my $info_gerp;	
	my $hit_position;
	my $GERP_SCORE;

	#I'm looking for rows with	
	#NCENC=TFM(VDR|VDR_JASPAR|chr7:5013558-5013573),TFM(VDR|VDR_XXmotif|chr7:5013558-5013573)
	#if there is no such line, there is no motifbr for vdr either
	next unless($_ =~ /TFM\(VDR\|VDR_JASPAR/); #TODO HARDCODED VDR_JASPAR
	my ($chr, $pos, $rs_id, $info) = (split /\t/)[0,1,2,7];
	#FILTER: BED
	if($subset_bed){ next unless(defined check_coords_in_bed($chr, $pos)); }
	$chr =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19
	my $genomic_coord = $chr . '-' . $pos;
	my @info = split(";", $info);
	#I want to flag rows which have
	#1 a NCENC field with a VDR_JASPAR in it
	#2 a MOTIFBR field with a VDR_JASPAR in it
	#3 both
	foreach my $item (@info){
		$info_ncenc = $item if( ($item =~ /^NCENC/)     &&  ($item =~ /TFM\(VDR\|VDR_JASPAR\|(.*)\)/) ); 
		$info_motifbr = $item if( ($item =~ /^MOTIFBR/) && ($item =~ /VDR_JASPAR/));
		$info_gerp = $item if($item =~ /GERP/);			
	}
	next unless($info_gerp);
	#FILTER: ENHANCER
	if($ENH_ONLY){ next unless($info_ncenc =~ /Enhancer/); }
	#get the GERP score
	if($info_gerp =~ /^GERP=(.*)/){
		$GERP_SCORE = $1;
		next if($GERP_SCORE eq '.'); #for some reason there are entries with no score and a dot
	}else{
		print STDERR "Warning: Gerp score: $info_gerp not recognised. Skipping.\n";
		next;
	}

	#looking for the hit position for this nucleotide, so that we can fill the hash
	if($info_motifbr){#===================================================================================================
		my ($motif_change_type,$motif_change_byTF) = split("=", $info_motifbr); 
		my @motif_change_byTF = split(",", $motif_change_byTF);
		
		foreach my $item (@motif_change_byTF){
			my ($TF,$TF_motif,$TF_motif_start,$TF_motif_end,$TF_motif_strand,$TF_snp_position,$TF_alt_score,$TF_ref_score) =  split("#", $item);
			next unless( $TF_motif eq $motif_name);
			
			next if(!$TF_motif_strand);
			next unless( ($TF_motif_strand eq '-') or ($TF_motif_strand eq '+') );
			
			$hit_position = get_hit_position($TF_motif_strand, $TF_snp_position);
			$motifmodel_motifpos2geomicpos2GERP_MOTIFBR{$hit_position}{$genomic_coord} = $GERP_SCORE; #at each genomic coordinate there can be only 1 gerp score
		}	
	}else{ #no motifbr - get hit position manually
		my $motif_coordinates;
		my $TF_motif_strand;
		my $TF_snp_position;
		
		my ($ncenc,$ncenc_data) = split("=", $info_ncenc);
		my @ncenc_data = split(",", $ncenc_data);
		
		foreach my $item (@ncenc_data){
			if($item =~ /TFM\(VDR\|VDR_JASPAR\|(.*)\)/){ #TODO HARDCODED VDR_JASPAR
				$motif_coordinates = $1;
				$TF_motif_strand = $vdr_motif_to_strand{$motif_coordinates};
				
				next if(!$TF_motif_strand);
				next unless( ($TF_motif_strand eq '-') or ($TF_motif_strand eq '+') );	
				
				#calculate snp position in motif
				my ($thischr,$coords) = split(":", $motif_coordinates); 
				my ($TF_motif_start, $TF_motif_end) = split('-', $coords);
				my $motif_size = $TF_motif_end - $TF_motif_start;# has to be 15
				if($motif_size ne $MOTIF_SIZE){
					print STDERR "Warning motif size: $motif_size is not 15. Skipping..\n";
					next;
				}
				$TF_snp_position =  $pos - $TF_motif_start; #pos and TF_motif start are absolute, this gives relative snp position
				$hit_position = get_hit_position($TF_motif_strand, $TF_snp_position);
				$motifmodel_motifpos2geomicpos2GERP{$hit_position}{$genomic_coord} = $GERP_SCORE;
			}
		}
	}	
}
close $instream;

#out format:
#type	position	score
#gob	pos_1	x
#lob	pos_1 y
#motifbr	pos_1 z
#...
my $TYPE_BR;my $TYPE_NONBR;
if($sym_variants){
#	$TYPE_BR = 'MOTIFBR_SYM';
#	$TYPE_LOB = 'LOB_SYM';
#	$TYPE_GOB = 'GOB_SYM';
	$TYPE_BR    = 'SYM';
	$TYPE_NONBR = 'SYM';
}else{
	$TYPE_BR    = 'MOTIFBR';
	$TYPE_NONBR = 'NO_MOTIFBR';
}

open (my $outstream,  q{>}, $output) or die("Unable to open $output : $!");
print $outstream "type\tposition\tGERP\n";
#motifbr
foreach my $motif_pos (@MOTIF_POSITIONS){
	if($motifmodel_motifpos2geomicpos2GERP_MOTIFBR{$motif_pos}){
		foreach my $genomic_pos (sort keys %{ $motifmodel_motifpos2geomicpos2GERP_MOTIFBR{$motif_pos} }){
			print $outstream "$TYPE_BR\t$motif_pos\t$motifmodel_motifpos2geomicpos2GERP_MOTIFBR{$motif_pos}{$genomic_pos}\n";
		}		
	}else{
		print $outstream "$TYPE_BR\t$motif_pos\tNA\n";
	}
}
#no motif br
foreach my $motif_pos (@MOTIF_POSITIONS){
	if($motifmodel_motifpos2geomicpos2GERP{$motif_pos}){
		foreach my $genomic_pos (sort keys %{ $motifmodel_motifpos2geomicpos2GERP{$motif_pos} }){
			print $outstream "$TYPE_NONBR\t$motif_pos\t$motifmodel_motifpos2geomicpos2GERP{$motif_pos}{$genomic_pos}\n";
		}		
	}else{
		print $outstream "$TYPE_NONBR\t$motif_pos\tNA\n";
	}
}
close $outstream;











#-------------------------------------subs-------------------------------
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

sub get_hit_position{
	my ($motif_strand, $snp_position) = @_;
	my $output;

	if($motif_strand =~ /-/){
		my $snp_position_RC = ($MOTIF_SIZE - $snp_position + 1 );
		$output = 'pos_' . $snp_position_RC;
	}elsif($motif_strand =~ /\+/){
		$output = 'pos_' . $snp_position;
	}else{
		print STDERR "get_hit_position() - motif hit strand: $motif_strand not recognised, skipping..\n";
		return undef;
	}
	return $output;
}