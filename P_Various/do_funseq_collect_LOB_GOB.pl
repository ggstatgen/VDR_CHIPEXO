#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#23/3/2015 Chris wants me to discard all positions that do not have ancestral state. Maybe count them to see how many are there?

#8/1/2015 This uses ALL VDR-BVS intersecting any VDR:RXR motif instance - not only MOTIFBRs!! Therefore for most of these you won't have a motif break score

#14/1/2015: Chris P: label those that break the motif in some way DONE
#14/1/2015: Make it possible to only select instances intersecting class I or II VDR peaks DONE 
#14/1/2015: make loss of binding negative and gain of binding positive: you want to show log_2 [ read_count (alt) / read_count(ref)] DONE
#14/1/2015: remove the box plot, leave only the median  - DO IN R

#17/1/2015:
#I will let R carry out the division.
#I will save the value (ANC)-(DER)

#Output is as follows:
#type	position	ref	alt

#I want to plot how the ASB-snps intersecting the VDRRXR Jaspar motif correlate with the phenotype - globally.
#I also want to see if the SYM-snps (those which are not causing significant asymmentric VDR binding) intersecting the VDRRX motif correlate with phenotype. 
#They should not correlate 
#The phenotype here is a measure of the difference between maternal and paternal reads as assessed by alleleseq

#For example, if you do the log2() of the ratio between ancestral and derived:
#+ samples will be cases where the derived provokes binding disruption (loss of binding variant - LOB)
#- samples will be cases where the derived allows binding (gain of binding variant - GOB)

#the current model is
#log2( (DER_count + 1) / (ANC_count + 1) ) < 0  => LOB
#log2( (DER_count + 1) / (ANC_count + 1) ) > 0  => GOB

#this is similar to do_funseq_build_motifdisruption_model.pl but it looks at the correlation of IMPACT of variant on VDR:RXR motif position WITH PHENOTYPE (binding affinity)

#Input:
#1)######################################################
#Output.vcf file from funseq - I need lines tagged with NCENC and TFM VDR_JASPAR (labelled as falling under the JASPAR motif)

#2)######################################################
#all raw alleleseq outputs with raw count reads from REF carrying allele parent  and ALT-carrying allele parent, for each sample
#typical input
#chrm	snppos   	ref	mat_gtyp	pat_gtyp	c_gtyp	phase	mat_all	pat_all	cA	cC	cG	cT	winning	SymCls	SymPval	BindingSite	cnv
#1	121428238	C	M	Y	M	PHASED	A	C	13	3	0	0	M	Asym	0.021270752	0	1
#1	121428255	T	M	Y	W	PHASED	A	T	19	0	0	0	M	Asym	3.81469726562E-006	0	1

#3)#########################################################
#ALL REFERENCES to REF ALLELE mean ANCESTRAL ALLELE to me (given my funseq analysis)
#So the third input is the full set of 1kg snps tested for ASB/QTL, from which I will get ancestral information

#4)####################################################
#The motif file I originally used for the funseq analysis. I need these to see if the motif is on the + or - strand, in order to compute WHERE
#in this motif the SNP is hitting (position)

#general eg
#score lnFC = ln [ (anc_count + 1) / (alt_count + 1) ]
#eg an asym site has 7 reads mapped to the alternate and 0 ancestral  => lnFC = ln[ 1 / 8  ]  = -3
# eg an asym site has 2 reads mapped to  the alternate and 9 to the ancestral => lnFC = ln[ 10 / 3 ] = 1.74

#for each position gather these mapped read values for all the 20 samples and build 20 lnFC values. Save these in histogram

#OUTPUT: #two files as follows:
#a matrix having a ROW per VDR:RXR position (15 column)
#one for LOB cases, one for GOB cases
#you WILL need to postprocess these two in gedit to order rows, replace separators with tabs, TRANSPOSE (in libreoffice calc)
#then you can plot in R, on the same graph, with different dots.

my $infile;
my $INPUT_ALLELESEQ_RAW = "/lustre/scratch110/sanger/gg14/VDR/ALLELESEQ/RESULTS/_RAW/interestingHets_vdrchipexo.txt";
my $INPUT_VARIANTS = '/lustre/scratch110/sanger/gg14/VDR/VARIANTS/ALLELESEQ_20SAMPLES/ALLELESEQ_20samples_merged.vcf.gz';
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
#$infile = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/Output.vcf";
#$infile = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_sym_allsamples_ancestral/Output.vcf";
#$ENH_ONLY = 1;
#$subset_bed = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_motif_scan_pscanchip/Pscanchip_occurrences_VDRRXR_jaspar_classI.bed";
#$sym_variants = 1;

if(!$infile){
	print "USAGE: do_funseq_collect_LOB_GOB.pl -i=<INFILE> -subset=<SUBSET_BED> -sym -enh\n";
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
	$output  = $directory . $basename . '_phdisruption_subsetbed.Rdata';
}elsif($ENH_ONLY){
	$output  = $directory . $basename . '_phdisruption_enhonly.Rdata';
}else{
	$output  = $directory . $basename . '_phdisruption_all.Rdata';
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
#1 -build huge hash to map chr-pos to ref, ancestral and alternate
#####################
my %variants_1kg;
my $counter_multiallelic = 0;
my $counter_na = 0;
my $counter_noancestral = 0;
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
	my $ref = $fields[3]; 
	my $alt = $fields[4];

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
#			$variants_1kg{$key}{ANC} = 'NA';
#			$variants_1kg{$key}{DER} = 'NA';
			$counter_na++;	
			next;
		}else{ 
			#the ancestral allele is a base, but a different one from ref and alt
			#so both ref and alt are "derived"
#			$variants_1kg{$key}{ANC} = 'NA';
#			$variants_1kg{$key}{DER} = 'NA';
			#Question asked to 1kg staff
			#Ignore these cases. The primates all had an allele and then a variation happened		
			$counter_multiallelic++;
			next;
		}
	}else{
#		$variants_1kg{$key}{ANC} = 'NA';
#		$variants_1kg{$key}{DER} = 'NA';	
		$counter_noancestral++;
		next
	}

#	#get ancestral
#	my @info = split(";", $fields[7]);
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
print STDERR "Multiallelic cases: $counter_multiallelic\n";
print STDERR "No ancestral cases: $counter_noancestral\n";
print STDERR "Cases with undefined ancestral: $counter_na\n";
#This should have taken a long time. now you know what the ancestral reference and alternate reference is.

#####################
#2 get all the motif instance coordinates with their orientation.
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


###################
#3 build interestinghets hash
###################
#chr-position => samplename => alleles / read counts info
my %position2sample2readdepth;
open ($instream,  q{<}, $INPUT_ALLELESEQ_RAW) or die("Unable to open $INPUT_ALLELESEQ_RAW : $!");
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
	
#	my $allele_ancestral = $ref;
#	my $allele_derived;
#	if($ref eq $m_allele){
#		$allele_derived = $p_allele;
#	}elsif($ref eq $p_allele){
#		$allele_derived = $m_allele;
#	}else{
#		next;
#	}

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
#4 pick alleleseq snps hitting motifs from funseq
#also pick motifBRs and label them as a third category
#########################
my %motifmodel_motifpos2geomicpos2FC_LOB;
my %motifmodel_motifpos2geomicpos2FC_GOB;
my %motifmodel_motifpos2geomicpos2FC_MOTIFBR;
open ($instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/);
	next if($_ eq '');
	#I'm looking for rows with	
	#NCENC=TFM(VDR|VDR_JASPAR|chr7:5013558-5013573),TFM(VDR|VDR_XXmotif|chr7:5013558-5013573)
	#if there is no such line, there is no motifbr for vdr either
	next unless($_ =~ /TFM\(VDR\|VDR_JASPAR/); #TODO HARDCODED VDR_JASPAR
	
	my $info_ncenc;
	my $info_motifbr;
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
	}
	#FILTER: ENHANCER
	if($ENH_ONLY){ next unless($info_ncenc =~ /Enhancer/); }
	my $hit_position;
	
	if($info_motifbr){#===================================================================================================
		my ($motif_change_type,$motif_change_byTF) = split("=", $info_motifbr); my @motif_change_byTF = split(",", $motif_change_byTF);
		
		foreach my $item (@motif_change_byTF){
			my ($TF,$TF_motif,$TF_motif_start,$TF_motif_end,$TF_motif_strand,$TF_snp_position,$TF_alt_score,$TF_ref_score) =  split("#", $item);
			next unless( $TF_motif eq $motif_name);
			next if(!$TF_motif_strand);
			next unless( ($TF_motif_strand eq '-') or ($TF_motif_strand eq '+') );
			
			if(!$position2sample2readdepth{$genomic_coord}){
				print STDERR "$rs_id - $genomic_coord: no ancestral/derived info available for this VDR-BV. Skipping..\n";
				next;
			}		
			
			$hit_position = get_hit_position($TF_motif_strand, $TF_snp_position);
			my @rd_fc_MOTIFBR; my $fc_MOTIFBR_this_asb_snp;
			#get read depth vector at this position for all samples
			foreach my $sample (sort keys %{ $position2sample2readdepth{$genomic_coord} }){
				my ($values, $fc) = process_read_data($position2sample2readdepth{$genomic_coord}{$sample}, $sample, $rs_id);
				if($fc){
					push(@rd_fc_MOTIFBR, $values);
					$fc_MOTIFBR_this_asb_snp = join(",",@rd_fc_MOTIFBR);
					$motifmodel_motifpos2geomicpos2FC_MOTIFBR{$hit_position}{$genomic_coord} = $fc_MOTIFBR_this_asb_snp;
				}else{
					print "Warning fold change value $fc is empty or undefined. Skipping..\n";
					next;
				}
			}
		}	
	}else{ #no motifbr - get hit and coordinates manually =================================================================
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
				if(!$position2sample2readdepth{$genomic_coord}){
					print STDERR "$rs_id - $genomic_coord: no ancestral/derived info available for this VDR-BV. Skipping..\n";
					next;
				}	
				
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
						
				my @rd_fc_LOB; my @rd_fc_GOB;
				my $fc_LOB_this_asb_snp; my $fc_GOB_this_asb_snp;
				#get read depth string for this ASB SNP from alleleseq and process it into two array of fc values (if >= 1 or < 1)
				foreach my $sample (sort keys %{ $position2sample2readdepth{$genomic_coord} }){
					my ($values, $fc) = process_read_data($position2sample2readdepth{$genomic_coord}{$sample}, $sample, $rs_id);
					if($fc >= 1){ #TODO adjust the case for = 
						push(@rd_fc_LOB, $values); #the alternate allele produces loss of binding compared to the ref
						$fc_LOB_this_asb_snp = join(",",@rd_fc_LOB);
						$motifmodel_motifpos2geomicpos2FC_LOB{$hit_position}{$genomic_coord} = $fc_LOB_this_asb_snp;
					}elsif($fc < 1){
						push(@rd_fc_GOB, $values); #the alternate allele increases binding compared to the ref
						$fc_GOB_this_asb_snp = join(",",@rd_fc_GOB);
						$motifmodel_motifpos2geomicpos2FC_GOB{$hit_position}{$genomic_coord} = $fc_GOB_this_asb_snp;
					}else{
						print "Warning fold change value $fc is empty or undefined. Skipping..\n";
						next;
					}
				}				
			}
		}
	}	
}
close $instream;

#out format:
#type	position	ref	alt
#gob	pos_1 1.45
#lob	pos_1 1.5
#motifbr	pos_1	1.1
#...

my $TYPE_BR;my $TYPE_LOB; my $TYPE_GOB;
if($sym_variants){
	$TYPE_BR  = 'MOTIFBR_SYM';
	$TYPE_LOB = 'LOB_SYM';
	$TYPE_GOB = 'GOB_SYM';
#	$TYPE_BR  = 'SYM';
#	$TYPE_LOB = 'SYM';
#	$TYPE_GOB = 'SYM';	
}else{
	$TYPE_BR  = 'MOTIFBR';
	$TYPE_LOB = 'LOB';
	$TYPE_GOB = 'GOB';
}

open (my $outstream,  q{>}, $output) or die("Unable to open $output : $!");
print $outstream "type\tposition\tref\talt\n";

#motifbr
foreach my $motif_pos (@MOTIF_POSITIONS){
	if($motifmodel_motifpos2geomicpos2FC_MOTIFBR{$motif_pos}){
		my $all_fc_values_for_this_motif_pos = '';
		foreach my $genomic_pos (sort keys %{ $motifmodel_motifpos2geomicpos2FC_MOTIFBR{$motif_pos} }){
			my @fc = split(",", $motifmodel_motifpos2geomicpos2FC_MOTIFBR{$motif_pos}{$genomic_pos});
			foreach my $fc_value (@fc){
				my ($val_ref, $val_alt) = split('-', $fc_value);
				#print $outstream "$TYPE_BR\t$motif_pos\t$fc_value\n";
				print $outstream "$TYPE_BR\t$motif_pos\t$val_ref\t$val_alt\n";
			}
		}		
	}else{
		print $outstream "$TYPE_BR\t$motif_pos\tNA\tNA\n";
	}
}
#lob
foreach my $motif_pos (@MOTIF_POSITIONS){
	if($motifmodel_motifpos2geomicpos2FC_LOB{$motif_pos}){
		my $all_fc_values_for_this_motif_pos = '';
		foreach my $genomic_pos (sort keys %{ $motifmodel_motifpos2geomicpos2FC_LOB{$motif_pos} }){
			my @fc = split(",", $motifmodel_motifpos2geomicpos2FC_LOB{$motif_pos}{$genomic_pos});
			foreach my $fc_value (@fc){
				my ($val_ref, $val_alt) = split('-', $fc_value);
				#print $outstream "$TYPE_LOB\t$motif_pos\t$fc_value\n";
				print $outstream "$TYPE_LOB\t$motif_pos\t$val_ref\t$val_alt\n";
			}
		}		
	}else{
		print $outstream "$TYPE_LOB\t$motif_pos\tNA\tNA\n";
	}
}
#gob
foreach my $motif_pos (@MOTIF_POSITIONS){
	if($motifmodel_motifpos2geomicpos2FC_GOB{$motif_pos}){
		my $all_fc_values_for_this_motif_pos = '';
		foreach my $genomic_pos (sort keys %{ $motifmodel_motifpos2geomicpos2FC_GOB{$motif_pos} }){
			my @fc = split(",", $motifmodel_motifpos2geomicpos2FC_GOB{$motif_pos}{$genomic_pos});
			foreach my $fc_value (@fc){
				my ($val_ref, $val_alt) = split('-', $fc_value);
				#print $outstream "$TYPE_GOB\t$motif_pos\t$fc_value\n";
				print $outstream "$TYPE_GOB\t$motif_pos\t$val_ref\t$val_alt\n";
			}
		}		
	}else{
		print $outstream "$TYPE_GOB\t$motif_pos\tNA\tNA\n";
	}
}
close $outstream;


#subroutine===============================================================================
#this splits a read_data string of the kind "c,a,,0,6,0,0"
#into a phenotype disruption score FC = [ (alt_count + 1) / (anc_count + 1) ] (turned into logs in R)
#I do alt/anc so that LOB points are negative on the graph, while GOB are positive
#if the motif hitting snp with the alternative allele has lowered the coverage at the position, this subroutine should output a negative value

#input meaning:
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
	
	#I will let R carry out the division.
	#I will save the value (REF)-(ALT)
	my $value_string = $anc_count . '-' . $der_count;
	return ($value_string, $fold_change);
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
