#!/usr/bin/perl
use strict;
use warnings;
use List::MoreUtils qw(uniq);
use Getopt::Long;
use File::Basename;
use IO::Zlib;

#This only collects DAF in the VDR:RXR Jaspar motif. I used it to produce the data in figure 4C. I only selected the data in positions 1-6 ad 10-15

#23/3/2015 Chris wants me to discard all positions that do not have ancestral state. Maybe count them to see how many are there?

#25/2/2015 
#I collected and plotted phastcons scores for each nucleotide at the motif positions.
#This gives me a perspective over vertebrate conservation at the motif

#I now want to see what happens at the population level: I will collect derived allele frequencies for the VDR-BV variants.
#The hypothesis is that these are different based on the nucleotide

#info from 1000genomes
##INFO=<ID=AC,Number=.,Type=Integer,Description="Alternate Allele Count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total Allele Count">
##INFO=<ID=AMR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from AMR based on AC/AN">
##INFO=<ID=ASN_AF,Number=1,Type=Float,Description="Allele Frequency for samples from ASN based on AC/AN">
##INFO=<ID=AFR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from AFR based on AC/AN">
##INFO=<ID=EUR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from EUR based on AC/AN">

#ATTENTION: when the ancestral allele is the alternate allele, the info above is the ANCESTRAL ALLELE COUNT.
#To get the DERIVED ALLELE COUNT, do a 1-ANCESTRAL ALLELE COUNT

#FEATURE REQUEST: TAG all DAF datapoints according to their LOB/GOB potential
#however, now I'm showing datapoints. LOB and GOB are phenotypes. There is one per sample, not one per position. 
#Plot only when, if multiple, both agree.
#That is, the data to save is a DAF, with a phenotype label attached: this can be LOB or GOB, but not both. If it has both, don't plot it.

#http://www.1000genomes.org/category/frequently-asked-questions/population
#CEU are EUR, YRI are AFR

#input: funseq file, variants tested 
#to increase power only compare variants in hexamers to variants in dr3 space

#algorithm
#-save 1kg variants in hash. Data to keep is AFR_AF and EUR_AF

my $infile;
my $INPUT_MOTIF_FILE = '/net/isi-scratch/giuseppe/tools/funseq2-1.0/data/ENCODE.tf.bound.union.bed'; #hg19
my $INPUT_VARIANTS = '/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/ALLELESEQ_20samples_merged.vcf.gz'; #b37
my $INPUT_ALLELESEQ_RAW = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/_RAW/interestingHets_vdrchipexo.txt"; #b37
my $motif_name = 'VDR_JASPAR';
my $ENH_ONLY;
my $subset_bed; #if you only want to collect data for positions in an external bed file (eg CLASS I sites)
my $sym_variants; #flag, set it if you're dealing with sym variants

my @MOTIF_POSITIONS = ("pos_1", "pos_2", "pos_3", "pos_4", "pos_5", "pos_6", "pos_7", "pos_8", "pos_9", "pos_10", "pos_11", "pos_12", "pos_13", "pos_14", "pos_15");
my $MOTIF_SIZE = 15;
GetOptions(
        'i=s'      =>\$infile,
        'subset=s' =>\$subset_bed,
        'sym'      =>\$sym_variants,
        'enh'      =>\$ENH_ONLY
); 
#$infile = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/Output.vcf"; #hg19
#$subset_bed = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_motif_scan_pscanchip/Pscanchip_occurrences_VDRRXR_jaspar_classI_peakintervals.bed"; #hg19

if(!$infile){
	print "USAGE: do_funseq_collect_DAF.pl -i=<INFILE> -subset=<SUBSET_BED> -sym -enh\n";
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
	$output  = $directory . $basename . '_DAF_subsetbed.Rdata';
}elsif($ENH_ONLY){
	$output  = $directory . $basename . '_DAF_enhonly.Rdata';
}else{
	$output  = $directory . $basename . '_DAF_all.Rdata';
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

#####################
#2 -build  hash to map chr-pos to DAF-CEU and DAF-YRI
#####################
#before saving the frequencies, you need to to KNOW if the ancestral allele is the ref or the alt
#if the ancestral allele is the ref, save the frequency as is
#if the ancestral allele is the alt, the frequency you have is for the ancestral. The derived will be (1 - freq)
my %variants_1kg;
tie *FILE,   'IO::Zlib', $INPUT_VARIANTS, "rb";
while (<FILE>)	{ 
	chomp;
	next if($_ eq '');
	next if($_ =~ /^#/);
	my $ALTF_AFR; my $ALTF_EUR; my $ALTF_ASN; my $ALTF_AMR; #ALTERNATE allele frequencies
	my $FLAG; # set to one if the ancestral is the alternate	

	my @fields = split("\t", $_);
	#skip indels
	next if(length($fields[3]) > 1);
	next if(length($fields[4]) > 1);
	my $key = 'chr' . $fields[0] . '-' . $fields[1];
	#get derived allele frequencies
	my $ref = $fields[3]; 
	my $alt = $fields[4];
	my $info = $fields[7];

	#get ancestral allele info=====================
	my @info = split(";", $fields[7]);
	if($info =~ /AA=(\w+)/){
		my $anc = uc($1);
		if($anc eq $alt){
			$FLAG = 1;
			$variants_1kg{$key}{ANC} = $alt;
			$variants_1kg{$key}{DER} = $ref;	
		}elsif($anc eq $ref){
			$variants_1kg{$key}{ANC} = $ref;
			$variants_1kg{$key}{DER} = $alt;			
		}else{ 
			next;		
		}
	}else{
		next;	
	}	
	#get allele frequencies========================
	$ALTF_AFR = $1 if($info =~ /AFR_AF=([0-9]+\.[0-9]+)/);
	$ALTF_EUR = $1 if($info =~ /EUR_AF=([0-9]+\.[0-9]+)/);	
	$ALTF_AMR = $1 if($info =~ /AMR_AF=([0-9]+\.[0-9]+)/);
	$ALTF_ASN = $1 if($info =~ /ASN_AF=([0-9]+\.[0-9]+)/);		
	
	if($ALTF_EUR){
		if($FLAG){
			$variants_1kg{$key}{DAFEUR} = (1 - $ALTF_EUR);
		}else{
			$variants_1kg{$key}{DAFEUR} = $ALTF_EUR;
		}
	}
	if($ALTF_AFR){
		if($FLAG){
			$variants_1kg{$key}{DAFAFR} = (1 - $ALTF_AFR);
		}else{
			$variants_1kg{$key}{DAFAFR} = $ALTF_AFR;
		}
	}
	if($ALTF_AMR){
		if($FLAG){
			$variants_1kg{$key}{DAFAMR} = (1 - $ALTF_AMR);
		}else{
			$variants_1kg{$key}{DAFAMR} = $ALTF_AMR;
		}
	}
	if($ALTF_ASN){
		if($FLAG){
			$variants_1kg{$key}{DAFASN} = (1 - $ALTF_ASN);
		}else{
			$variants_1kg{$key}{DAFASN} = $ALTF_ASN;
		}
	}
}
close FILE;

###################
#3 get affinity binding phenotype info 
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
	my $coordinate_id = 'chr' . $chr . '-' . $snppos;
	
	my $ancestral = $variants_1kg{$coordinate_id}{ANC}; 
	my $derived = $variants_1kg{$coordinate_id}{DER}; 
	next unless($ancestral);
	next unless($derived);

	if( ($m_allele eq 'None') or ($p_allele eq 'None')  ){
	}else{
		if( ($ancestral eq $m_allele)  && ($derived eq $p_allele) ){		
		}elsif( ($ancestral eq $p_allele)  && ($derived eq $m_allele)  ){	
		}else{
			print STDERR "$coordinate_id: (ancestral ref/alt and mat/pat don't match: ($ancestral, $derived) and ($m_allele, $p_allele)\n";
			next;
		}		
	}
	my $read_data = join(",", $ancestral,$derived, $cA, $cC, $cG, $cT);
	$position2sample2readdepth{$coordinate_id}{$sample_id} = $read_data;
}
close $instream;

#for each motif position, for each coordinate, there might be multiple phenotypes.
#you want to attach one phenotype label per key
my %motifmodel_motifpos2geomicpos2DAF;
open ($instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/);
	next if($_ eq '');
	my $motif_coordinates;
	my $TF_motif_strand;
	my $TF_snp_position;
	#I'm looking for rows with	
	#NCENC=TFM(VDR|VDR_JASPAR|chr7:5013558-5013573),TFM(VDR|VDR_XXmotif|chr7:5013558-5013573)
	#if there is no such line, there is no motifbr for vdr either
	next unless($_ =~ /TFM\(VDR\|VDR_JASPAR/); #TODO HARDCODED VDR_JASPAR
	
	my $info_ncenc;	my $info_motifbr;
	my ($chr, $pos, $rs_id, $info) = (split /\t/)[0,1,2,7]; #chr is hg19
	#FILTER: BED
	if($subset_bed){ next unless(defined check_coords_in_bed($chr, $pos)); }
	
	#$chr =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19
	my $genomic_coord = $chr . '-' . $pos;
	my @info = split(";", $info);
	
	#I want to flag rows which have
	#1 a NCENC field with a VDR_JASPAR in it
	foreach my $item (@info){
		$info_ncenc = $item if( ($item =~ /^NCENC/)     &&  ($item =~ /TFM\(VDR\|VDR_JASPAR\|(.*)\)/) ); 
		$info_motifbr = $item if( ($item =~ /^MOTIFBR/) && ($item =~ /VDR_JASPAR/));		
	}
	#FILTER: ENHANCER
	if($ENH_ONLY){ next unless($info_ncenc =~ /Enhancer/); }
	my $hit_position;
	
	#if you have motifbr data, which you do, you can remove instances where motif pwm change is discordant with phenotype change.
	#if you do not have motifbr data, you cannot do this, and you will have to look at the final results
	
	#MOTIFBR=VDR#VDR_JASPAR#138803009#138803024#-#2#0.000#1.000
	if($info_motifbr){#===================================================================================================
		my ($motif_change_type,$motif_change_byTF) = split("=", $info_motifbr); my @motif_change_byTF = split(",", $motif_change_byTF);
		
		foreach my $item (@motif_change_byTF){
			my ($TF,$TF_motif,$TF_motif_start,$TF_motif_end,$TF_motif_strand,$TF_snp_position,$TF_der_score,$TF_anc_score) =  split("#", $item);
			next unless( $TF_motif eq $motif_name);
			next if(!$TF_motif_strand);
			next unless( ($TF_motif_strand eq '-') or ($TF_motif_strand eq '+') );
			
			if(!$position2sample2readdepth{$genomic_coord}){
				print STDERR "$rs_id - $genomic_coord: no ancestral/derived info available for this VDR-BV. Skipping..\n";
				next;
			}				
			
			my %concordant_samples; my %discordant_samples;
			$concordant_samples{$genomic_coord} = 0; $discordant_samples{$genomic_coord} = 0;
			#Compute gt_score from $TF_alt_score and $TF_alt_score
			#these should already account for ancestral allele in funseq
			my $GT_SCORE = ( $TF_anc_score + 1 ) / ($TF_der_score + 1); #should these all be  LOB by definition? And therefore GT_SCORE > 1 ?
			#Compute ph score and evaluate concordance with gt score at this position
			#if no concordance, it's meaningless to get the DAF
			foreach my $sample (sort keys %{ $position2sample2readdepth{$genomic_coord} }){
				my ($sample_values, $sample_fc) = process_read_data($position2sample2readdepth{$genomic_coord}{$sample}, $sample, $rs_id);
				if( ($GT_SCORE >= 1) && ($sample_fc >= 1) ){
					$concordant_samples{$genomic_coord} += 1;
				}elsif( ($GT_SCORE >= 1) && ($sample_fc < 1)  ){
					$discordant_samples{$genomic_coord} += 1;
				}elsif( ($GT_SCORE < 1) && ($sample_fc >= 1)  ){
					$discordant_samples{$genomic_coord} += 1;
				}elsif( ($GT_SCORE < 1) && ($sample_fc < 1)  ){
					print STDERR "$rs_id: motifbr gt-ph concordance, but GT_SCORE is negative?? Check.\n";
					$concordant_samples{$genomic_coord} += 1;
				}else{
					print STDERR "Should never be here\n";
				}
			}
			#Now I know, for each samples how many are concordant and how many are discordant. All are LOB, unless something strange is going on, but some will be false.
			if($discordant_samples{$genomic_coord} > 0){
				if($concordant_samples{$genomic_coord} > 0){
						print STDERR "Warning - motifbr for VDR-BV $rs_id shows:\n"; 
						print STDERR "concordance for $concordant_samples{$genomic_coord} samples and\n";
						print STDERR "discordance for $discordant_samples{$genomic_coord} samples.\n";
						print STDERR "Skipping (EVALUATE THIS)\n";
						next;
					next;
				}else{
					print STDERR "Warning - motifbr for VDR-BV $rs_id shows only phenotypes discordant ($discordant_samples{$genomic_coord}) with the genotype. Skipping..\n";
					next;
				}
			}
			#Save DAF info for this motif position and this coordinate
			$hit_position = get_hit_position($TF_motif_strand, $TF_snp_position);
			if($variants_1kg{$genomic_coord}{DAFEUR}){
				$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFEUR} = $variants_1kg{$genomic_coord}{DAFEUR};
			}else{
				$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFEUR} = 'NA';
			}
			if($variants_1kg{$genomic_coord}{DAFAFR}){
				$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFAFR} = $variants_1kg{$genomic_coord}{DAFAFR};
			}else{
				$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFAFR} = 'NA';
			}
			if($variants_1kg{$genomic_coord}{DAFAMR}){
				$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFAMR} = $variants_1kg{$genomic_coord}{DAFAMR};
			}else{
				$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFAMR} = 'NA';
			}			
			if($variants_1kg{$genomic_coord}{DAFASN}){
				$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFASN} = $variants_1kg{$genomic_coord}{DAFASN};
			}else{
				$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFASN} = 'NA';
			}
			#save LOB-MOTIFBR label
			$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{LABEL} = 'LOB,MOTIFBR';
		}	
	}else{ #no motifbr - get hit and coordinates manually 
		my $motif_coordinates;my $TF_motif_strand;my $TF_snp_position;
		my $LABEL;
		
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
				
				#get lob and gob label for all samples
				my %lob_samples; my %gob_samples; 
				$lob_samples{$genomic_coord} = 0; $gob_samples{$genomic_coord} = 0;
				foreach my $sample (sort keys %{ $position2sample2readdepth{$genomic_coord} }){
					my ($values, $fc) = process_read_data($position2sample2readdepth{$genomic_coord}{$sample}, $sample, $rs_id);
					if($fc > 1){
						$lob_samples{$genomic_coord} += 1;
					}elsif($fc < 1){
						$gob_samples{$genomic_coord} += 1;
					}elsif($fc == 1){
						print STDERR "Warning: VDR-BV $rs_id, sample $sample - phenotype fc is $fc. Skipping sample..\n";
						next;
					}else{
						print "Warning fold change value $fc is empty or undefined. Skipping..\n";
						next;
					}
				}
				
				#GET LOB/GOB labels
				if(  ($lob_samples{$genomic_coord} > 0 ) && ($gob_samples{$genomic_coord} > 0)  ){
					print STDERR "Warning VDR-BV $rs_id - at this position, $lob_samples{$genomic_coord} sample are LOB and $gob_samples{$genomic_coord} are GOB.\n";
					print STDERR "What to do? Skip for now.\n";
					next;
				}elsif( ($lob_samples{$genomic_coord} > 0 ) && ($gob_samples{$genomic_coord} == 0) ){
					$LABEL = 'LOB';
				}elsif( ($gob_samples{$genomic_coord} > 0 ) && ($lob_samples{$genomic_coord} == 0)  ){
					$LABEL = 'GOB';
				}else{
					print STDERR "ERROR VDR-BV $rs_id $lob_samples{$genomic_coord} lob samples and  $gob_samples{$genomic_coord} gob samples - should not be here.\n";
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
				
				#Save DAF info for this motif position and this coordinate
				$hit_position = get_hit_position($TF_motif_strand, $TF_snp_position);
				
				if($variants_1kg{$genomic_coord}{DAFEUR}){
					$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFEUR} = $variants_1kg{$genomic_coord}{DAFEUR};
				}else{
					$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFEUR} = 'NA';
				}
				if($variants_1kg{$genomic_coord}{DAFAFR}){
					$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFAFR} = $variants_1kg{$genomic_coord}{DAFAFR};
				}else{
					$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFAFR} = 'NA';
				}
				if($variants_1kg{$genomic_coord}{DAFAMR}){
					$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFAMR} = $variants_1kg{$genomic_coord}{DAFAMR};
				}else{
					$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFAMR} = 'NA';
				}			
				if($variants_1kg{$genomic_coord}{DAFASN}){
					$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFASN} = $variants_1kg{$genomic_coord}{DAFASN};
				}else{
					$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{DAFASN} = 'NA';
				}
				#save label
				$motifmodel_motifpos2geomicpos2DAF{$hit_position}{$genomic_coord}{LABEL} = $LABEL;
			}
		}
	}
}
close $instream;

#motifpos -> genomic_coord -> DAFCEU,...DAFASN
#motifpos -> genomic_coord -> LABEL

my $TYPE;
if($sym_variants){
	$TYPE = 'sym';
}else{
	$TYPE = 'asym';
}

open (my $outstream,  q{>}, $output) or die("Unable to open $output : $!");
print $outstream "TYPE\tPOSITION\tDIRECTION\tETHNICITY\tDAF\n";

foreach my $motif_pos (@MOTIF_POSITIONS){
	if($motifmodel_motifpos2geomicpos2DAF{$motif_pos}){
		foreach my $genomic_pos (sort keys %{ $motifmodel_motifpos2geomicpos2DAF{$motif_pos} }){
			my $direction = $motifmodel_motifpos2geomicpos2DAF{$motif_pos}{$genomic_pos}{LABEL};
			my $eur = $motifmodel_motifpos2geomicpos2DAF{$motif_pos}{$genomic_pos}{DAFEUR};
			my $afr = $motifmodel_motifpos2geomicpos2DAF{$motif_pos}{$genomic_pos}{DAFAFR};
			my $amr = $motifmodel_motifpos2geomicpos2DAF{$motif_pos}{$genomic_pos}{DAFAMR};
			my $asn = $motifmodel_motifpos2geomicpos2DAF{$motif_pos}{$genomic_pos}{DAFASN};
			print $outstream $TYPE . "\t" . $motif_pos . "\t" . $direction . "\t" . 'CEU' . "\t" . $eur . "\n";
			print $outstream $TYPE . "\t" . $motif_pos . "\t" . $direction . "\t" . 'YRI' . "\t" . $afr . "\n";
			print $outstream $TYPE . "\t" . $motif_pos . "\t" . $direction . "\t" . 'AMR' . "\t" . $amr . "\n";
			print $outstream $TYPE . "\t" . $motif_pos . "\t" . $direction . "\t" . 'ASN' . "\t" . $asn . "\n";
		}
	}else{
		print $outstream "$TYPE\t$motif_pos\tNA\tNA\tNA\n";
	}
}

close $outstream;

#============================================
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