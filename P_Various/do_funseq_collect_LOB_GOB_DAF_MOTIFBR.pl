#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#I have found interesting stuff about the DAF in position 12 (maybe not true anymore)
#Now I want to collect everything I know from:
#1 the funseq file
#2 the 1kg file
#3 the alleleseq file
#about phenotype variation (lob/gob), motifbreak yes/no, funseq annotation

#ONLY for the variants in one position specified

#INPUTS
#1kg vcf
#funseq output file
#alleleseq output file
#motif file

#my $INPUT_MOTIF_FILE = '/net/isi-scratch/giuseppe/tools/funseq2-1.0/data/ENCODE.tf.bound.union.bed'; #hg19
my $INPUT_MOTIF_FILE = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_motif_scan_pscanchip/Pscanchip_occurrences_VDRRXR.all.ris";
#I want the motif sequence and the motif score, or at least the score
my $INPUT_VARIANTS = '/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/ALLELESEQ_20samples_merged.vcf.gz'; #b37
#my $INPUT_VARIANTS = "/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/test_rs71404233.gz";
my $INPUT_ALLELESEQ_RAW = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/_RAW/interestingHets_vdrchipexo.txt"; #b37

#OUTPUT
#tsv, 1 row per variant with info containing: lob, gob, summary of funseq, daf
#TODO split funseq annotation and keep: GERP, HUB, NENC, HOT, GENE, RECUR
#TODO also plot pscanchip motif score

my $infile;
my $POSITION;
my $subset_bed;
my $sym_variants;			
my $ENH_ONLY;
my $MOTIF_SIZE = 15;
my $motif_name = 'VDR_JASPAR';
GetOptions(
        'i=s'      =>\$infile,
        'p=i'      =>\$POSITION,
        'subset=s' =>\$subset_bed,
        'sym'      =>\$sym_variants,
        'enh'      =>\$ENH_ONLY
);
#$infile = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/Output.vcf"; #hg19
#$POSITION = 12;
#$subset_bed = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_motif_scan_pscanchip/Pscanchip_occurrences_VDRRXR_jaspar_classI_peakintervals.bed"; #hg19

if(!$infile){
	print "USAGE: do_funseq_collect_LOB_GOB_DAF_MOTIFBR_bypos.pl -i=<INFILE> -p=<POSITION> -subset=<SUBSET_BED> -sym -enh\n";
    print "<INFILE> vcf output of funseq2\n";
    print "<POSITION> a number in the [1,15] range, corresponding with the nucleotide position of interest within VDR:RXR\n";
    print "(optional)<SUBSET_BED> only collect data for VDR-BVs breaking motifs intersecting with this bed file(default=no)\n";
    print "(optional)<sym> flag; set it if you are dealing with SYM variants from alleleseq (default=no)";
    print "(optional)<enh> whether to only use snps with NCENC=enhancer (default=no)\n";
    exit 1;
}
if(!$POSITION){
	print "USAGE: do_funseq_collect_LOB_GOB_DAF_MOTIFBR_bypos.pl -i=<INFILE> -p=<POSITION> -subset=<SUBSET_BED> -sym -enh\n";
    print "<INFILE> vcf output of funseq2\n";
    print "<POSITION> a number in the [1,15] range, corresponding with the nucleotide position of interest within VDR:RXR\n";
    print "(optional)<SUBSET_BED> only collect data for VDR-BVs breaking motifs intersecting with this bed file(default=no)\n";
    print "(optional)<sym> flag; set it if you are dealing with SYM variants from alleleseq (default=no)";
    print "(optional)<enh> whether to only use snps with NCENC=enhancer (default=no)\n";
    exit 1;
}
my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
my $output;

if($subset_bed){
	$output  = $directory . $basename . '_summary_LOB_GOB_DAF_MOTIFBR_pos_' . $POSITION . '_subsetbed.tsv';
}elsif($ENH_ONLY){
	$output  = $directory . $basename . '_summary_LOB_GOB_DAF_MOTIFBR_pos_' . $POSITION . '_enhonly.tsv';
}else{
	$output  = $directory . $basename . '_summary_LOB_GOB_DAF_MOTIFBR_pos_' . $POSITION . '_all.tsv';
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
#This is needed to find the position where the SNP hits in the motif sequence, the score and the motif sequence
my %vdr_motif_to_strand;
my %vdr_motif_to_score;
my %vdr_motif_to_sequence_plus; #I will reverse comp motifs on the - for compatibility with hit_position
open (my $instream,  q{<}, $INPUT_MOTIF_FILE) or die("Unable to open $INPUT_MOTIF_FILE : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^CHR/);
	my $motif_sequence;
	my @fields = split("\t", $_);
	#format is as follows
	#CHR     REG_START       REG_END REG_STRAND      ABS_SITE_START  ABS_SITE_END    REL_SITE_START  REL_SITE_END    SITE_STRAND     SCORE   SITE
	#chr9    131645230       131645379       +       131645239       131645253       -66     -53     +       0.989252        GGGTCATGGAGTTCA
	#chr6    34302932        34303081        +       34302996        34303010        -11     2       -       0.989252        TGAACTCCATGACCC
	
	#ris coordinates are -1 of the ones needed by UCSC. 
	#If you want the UCSC BROWSER ones, add 1 to both
	#If you want a standard BED, add 1 to end
	#If you want the BIGWIGTOWIG, add 1 to end
	
	#SAVING BIGWIGTOWIG coordinates
	my $start = $fields[4]; #bw2w will count from this + 1
	my $stop = $fields[5] + 1;
	my $identifier = $fields[0] . ':' . $start . '-' . $stop;    #chr7:5013558-5013573
	$vdr_motif_to_strand{$identifier} = $fields[8];
	$vdr_motif_to_score{$identifier} = $fields[9];
	
	$motif_sequence = $fields[10];
	if($fields[8] eq '-'){
		my $complement  = $motif_sequence;
		$complement =~ tr/ACGTacgt/TGCAtgca/;
		my $rc = reverse($complement);
		$motif_sequence = $rc;
	}
	$vdr_motif_to_sequence_plus{$identifier} = $motif_sequence;	
}
close $instream;



#####################
#2 -build hash to map chr-pos to ref, ancestral, alternate, and DAFs
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
	my $ref = uc($fields[3]); 
	my $alt = uc($fields[4]);
	my $info = $fields[7];
	
	#get ancestral allele info=====================
	my @info = split(";", $fields[7]);
	if($info[0] =~ /^AA=(.*)/){
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

	#reverse frequency if alternative frequency is ancestral frequency; save
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
	
#	if($variants_1kg{$key}{INFO}){
#		print STDERR "ATTENTION: positional overlap in 1kg hash: $key. What to do?\n";
#	}
#	#put anc and der as key?
	
	$variants_1kg{$key}{INFO} = $info;
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
	
#	my $ancestral_ref = $ref;
#	my $ancestral_alt;
#	if($ref eq $m_allele){
#		$ancestral_alt = $p_allele;
#	}elsif($ref eq $p_allele){
#		$ancestral_alt = $m_allele;
#	}else{
#		next;
#	}

	#if($ref ne $ancestral_ref){
	#	print STDERR "$coordinate_id: ref and ancestral differ. ref: $ref, ancestral_ref: $ancestral_ref, ancestral_alt: $ancestral_alt\n";
	#}
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



#get the variant's position in the motif. If not the position required, skip. If the position is the required one, print in output
my %data_entries;
open ($instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/);
	next if($_ eq '');
	#I'm looking for rows with	
	#NCENC=TFM(VDR|VDR_JASPAR|chr7:5013558-5013573)
	#if there is no such line, there is no motifbr for vdr either
	next unless($_ =~ /TFM\(VDR\|VDR_JASPAR/); #TODO HARDCODED VDR_JASPAR
	my $info_ncenc; my $info_motifbr;
	
	my ($chr, $pos, $rs_id, $info) = (split /\t/)[0,1,2,7];
	#FILTER: BED
	if($subset_bed){ next unless(defined check_coords_in_bed($chr, $pos)); }
	my $genomic_coord = $chr . '-' . $pos;
	if(!$position2sample2readdepth{$genomic_coord}){
		print STDERR "$rs_id - $genomic_coord: no ancestral/derived info available for this VDR-BV. Skipping..\n";
		next;
	}		
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
	my $hit_position; my $field_motif_br;
	if($info_motifbr){
		$field_motif_br = 'Y';
	}else{
		$field_motif_br = 'N';
	}
	my $motif_coordinates;my $TF_motif_strand;my $TF_snp_position;
	my ($ncenc,$ncenc_data) = split("=", $info_ncenc);
	my @ncenc_data = split(",", $ncenc_data);
	my $field_ph_dir;
	
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
			#check required motif position
			#next if($hit_position ne $POSITION);

			foreach my $sample (sort keys %{ $position2sample2readdepth{$genomic_coord} }){
				my ($field_ref_alt, $field_fc, $field_ref, $field_alt) = process_read_data($position2sample2readdepth{$genomic_coord}{$sample}, $sample, $rs_id);
				if($field_fc >= 1){ #LOB
					$field_ph_dir = 'LOB';
				}elsif($field_fc < 1){ #GOB
					$field_ph_dir = 'GOB';
				}else{
					$field_fc = '-';
					$field_ph_dir = '-';
				}
				
				if($sample =~ /^NA19(.*)/){
					$sample = $sample . '(YRI)';
				}else{
					$sample = $sample . '(CEU)';
				}
				
				#process the funseq info string
				my $field_funseq_anno = process_funseq_annotation($info);
				my $TF_motif_score = $vdr_motif_to_score{$motif_coordinates};
				my $TF_motif_seq_plus = $vdr_motif_to_sequence_plus{$motif_coordinates};
				my $line =  $chr       . "\t" . 
								 $pos       . "\t" . 
								 $rs_id     . "\t" . 
								 $field_ref . "\t" . 
								 $field_alt . "\t" . 
								 $sample    . "\t" . 
								 $field_motif_br  . "\t" . 
								 $TF_motif_strand . "\t" .
								 $TF_motif_score  . "\t" .
								 $TF_motif_seq_plus . "\t" .
								 $hit_position    . "\t" .
								 $field_ref_alt   . "\t" .
								 $field_ph_dir    . "\t" .
								 $field_fc        . "\t" . 
								 $variants_1kg{$genomic_coord}{DAFEUR}   . "\t" . 
								 $variants_1kg{$genomic_coord}{DAFAFR}   . "\t" .

							 	 $variants_1kg{$genomic_coord}{INFO} . "\t" . 									 
								 $field_funseq_anno;
			$data_entries{$line} = 1; 
			}				
		}
	}	
}
close $instream;


open (my $outstream,  q{>}, $output) or die("Unable to open $output : $!");
#print $outstream "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
print $outstream "CHR\tPOS\tID\tREF\tALT\tSAMPLE\tMOTIF_BREAK\tMOTIF_STRAND\tMOTIF_SCORE\tMOTIF_SEQ(+)\tMOTIF_HIT_POS\tRD_REF_ALT\tPH_DIR\tPH_FC\tDAF_EUR\tDAF_YRI\t1KG_ANNOTATION\tFUNSEQ_GERP\tFUNSEQ_PROTEIN_NET_HUB\tFUNSEQ_NON_CODING_ANNOTATION\tFUNSEQ_HOT_REG\tFUNSEQ_GENE\tFUNSEQ_RECUR\n";
foreach my $item (keys %data_entries){
	print $outstream $item, "\n";
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
	
	my ($ancestral_ref,$ancestral_alt,$cA,$cC,$cG,$cT) = split(",", $read_data_string);
	my $ref_count; my $alt_count;
	my $alt_base = $ancestral_alt;
	my $ref_base = $ancestral_ref;
	
	if( ($cA eq 0) && ($cC eq 0) && ($cG eq 0) && ($cT eq 0) ){
		print STDERR "process_read_data(): error asb snp: $id. (LCL sample: $lcl_sample) - 0 read coverage for all possible nucleotides at this event ($cA,$cC,$cG,$cT). Aborting..\n";
		exit -1;
	}

	#get reference read coverage
	if($ref_base eq 'A'){
		$ref_count = $cA;
	}elsif($ref_base eq 'C'){
		$ref_count = $cC;
	}elsif($ref_base eq 'G'){
		$ref_count = $cG;
	}elsif($ref_base eq 'T'){
		$ref_count = $cT;
	}else{
		print STDERR "process_read_data(): error asb snp: $id. (LCL sample: $lcl_sample) - reference base not recognised: $ref_base. Aborting..\n";
		exit -1;
	}
	
	if($alt_base eq 'A'){
		$alt_count = $cA;
	}elsif($alt_base eq 'C'){
		$alt_count = $cC;
	}elsif($alt_base eq 'G'){
		$alt_count = $cG;
	}elsif($alt_base eq 'T'){
		$alt_count = $cT;
	}else{
		print STDERR "process_read_data(): error asb snp: $id. (LCL sample: $lcl_sample) - alternate base not recognised: $alt_base. Aborting..\n";
		exit -1;
	}
	my $fold_change = ( $ref_count + 1 ) / ( $alt_count + 1 );
	
	#I will let R carry out the division.
	#I will save the value (REF)-(ALT)
	my $value_string = $ref_count . '-' . $alt_count;
	return ($value_string, $fold_change, $ref_base, $alt_base);
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
		$output = $snp_position_RC;
	}elsif($motif_strand =~ /\+/){
		$output = $snp_position;
	}else{
		print STDERR "get_hit_position() - motif hit strand: $motif_strand not recognised, skipping..\n";
		return undef;
	}
	return $output;
}

#I have something like
#SAMPLE=interestingHets_NA19213_EBLfiltered_hg19;GERP=-3.66;CDS=No;HUB=FAM200B:REG(0.880),FBXL5:PPI(0.361);NCENC=DHS(MCV-42|chr4:15755985-15756135),Enhancer(chmm/segway|chr4:15755162-15756284),Enhancer(drm|chr4:15756000-15756200),TFM(VDR|VDR_JASPAR|chr4:15756098-15756113),TFM(VDR|VDR_XXmotif|chr4:15756098-15756113),TFM(VDR|VDR_dreme|chr4:15756098-15756105),TFP(FOS|chr4:15755301-15756613),TFP(GATA2|chr4:15755241-15756612),TFP(MAX|chr4:15755846-15756384),TFP(MYC|chr4:15755856-15756339),TFP(VDR|chr4:15755910-15756250);MOTIFBR=VDR#VDR_dreme#15756098#15756105#-#4#0.000000#0.476400,VDR#VDR_JASPAR#15756098#15756113#-#12#0.000#1.000,VDR#VDR_XXmotif#15756098#15756113#-#12#0.02457#0.73562;GENE=FAM200B(Distal)[pearson(H3K27ac):0.918602,pearson(H3K4me1):0.912105],FBXL5(Distal)[pearson(H3K27ac):0.875722,pearson(H3K4me1):0.873539];NCDS=1.67156771704517

#or

#SAMPLE=interestingHets_NA19248_EBLfiltered_hg19;GERP=0.375;CDS=No;NCENC=DHS(MCV-31|chr7:7298040-7298190),Enhancer(chmm/segway|chr7:7296800-7298400),TFM(VDR|VDR_JASPAR|chr7:7298122-7298137),TFP(BATF|chr7:7297738-7298801),TFP(BCL11A|chr7:7297625-7298314),TFP(IRF4|chr7:7297670-7298802),TFP(MAX|chr7:7297536-7298947),TFP(MEF2A|chr7:7297922-7298725),TFP(MYC|chr7:7297025-7298988),TFP(NFKB1|chr7:7297551-7298769),TFP(PAX5|chr7:7297634-7298780),TFP(PAX5|chr7:7297644-7298634),TFP(PAX5|chr7:7297659-7298342),TFP(PAX5|chr7:7297918-7298625),TFP(POU2F2|chr7:7297992-7298671),TFP(RXRA|chr7:7297871-7298750),TFP(SPI1|chr7:7297558-7298322),TFP(SPI1|chr7:7297604-7298306),TFP(VDR|chr7:7298055-7298256);HOT=Gm12878;MOTIFBR=VDR#VDR_JASPAR#7298122#7298137#+#12#0.000#1.000;NCDS=1.788999949

#SAMPLE=interestingHets_NA19213_EBLfiltered_hg19;GERP=-3.66;CDS=No;HUB=FAM200B:REG(0.880),FBXL5:PPI(0.361);NCENC=DHS(MCV-42|chr4:15755985-15756135),Enhancer(chmm/segway|chr4:15755162-15756284),Enhancer(drm|chr4:15756000-15756200),TFM(VDR|VDR_JASPAR|chr4:15756098-15756113),TFM(VDR|VDR_XXmotif|chr4:15756098-15756113),TFM(VDR|VDR_dreme|chr4:15756098-15756105),TFP(FOS|chr4:15755301-15756613),TFP(GATA2|chr4:15755241-15756612),TFP(MAX|chr4:15755846-15756384),TFP(MYC|chr4:15755856-15756339),TFP(VDR|chr4:15755910-15756250);MOTIFBR=VDR#VDR_dreme#15756098#15756105#-#4#0.000000#0.476400,VDR#VDR_JASPAR#15756098#15756113#-#12#0.000#1.000,VDR#VDR_XXmotif#15756098#15756113#-#12#0.02457#0.73562;GENE=FAM200B(Distal)[pearson(H3K27ac):0.918602,pearson(H3K4me1):0.912105],FBXL5(Distal)[pearson(H3K27ac):0.875722,pearson(H3K4me1):0.873539];NCDS=1.67156771704517


#I want
#FUNSEQ_GERP	FUNSEQ_PROTEIN_NET_HUB	FUNSEQ_NON_CODING ANNOTATION	FUNSEQ_HOT_REG	FUNSEQ_GENE	FUNSEQ_RECUR
#GERP=-1.03	HUB=HTT:PPI(0.902)	DHS, Enhancer, ENCODE TFP	Y	HTT intron	
#GERP=-6.13	HUB=FOXO1:PPI(0.787)	DHS, Enhancer, ENCODE TFP	Y	FOXO1 intron	
#GERP=2.51	HUB=TLE3:PHOS(0.459)PPI(0.448)REG(0.761)	Enhancer, ENCODE TFP	Y	NOX5 distal, TLE3 distal	NA19190 (YRI)
#GERP=-3.66	HUB=FAM200B:REG(0.880),FBXL5:PPI(0.361)	DHS, Enhancer, ENCODE TFP		FAM200B (Distal), FBXL5(Distal)	
#GERP=0.375	NA	DHS, Enhancer, ENCODE TFP	Y		
#GERP=3.39	HUB=SCAND1:PPI(0.657)REG(0.959)	ENCODE TFP	Y	SCAND1 (Medial&UTR)	


sub process_funseq_annotation{
	my ($funseq_string) = @_;
	#output fields
	my $info_gerp; my $info_hub; my $info_ncenc; my $info_hot; my $info_gene; my $info_recur;
	my $out_gerp; my $out_hub; my $out_ncenc; my $out_hot; my $out_gene; my $out_recur;
	my @info = split(";",$funseq_string);
	
	foreach my $item (@info){
		$info_gerp = $item if( ($item =~ /^GERP/) );
		$info_hub = $item if( ($item =~ /^HUB/) ) ;				
		$info_ncenc = $item if( ($item =~ /^NCENC/) );
		$info_hot = $item if( ($item =~ /^HOT/) );
		$info_gene = $item if( ($item =~ /^GENE/) );
		$info_recur = $item if( ($item =~ /^RECUR/) );		
	}
	#gerp-----------------------------------------------------
	if($info_gerp =~ /GERP=(.*)/){
		$out_gerp = $1;
	}else{
		$out_gerp = '';
	}
	
	#hub-----------------------------------------------------
	if($info_hub){
		my %entries; my @entries;
		my ($header, $data) = split('=',$info_hub);
		my @info_hub = split(',', $data);
		foreach my $item (@info_hub){
			if($item =~ /(.+)\:(.+)/){
				$entries{$1} = 1;
			}else{
				next;
			}				
		}
		foreach my $item (sort keys %entries){
			push(@entries, $item);	
		}
		$out_hub = join(',',@entries);
	}else{
		$out_hub = '';
	}
	
	#ncenc-----------------------------------------------------
	if($info_ncenc){
		my %entries; my @entries;
		my ($header, $data) = split('=',$info_ncenc);
		my @info_ncenc = split(',', $data);
		foreach my $item (@info_ncenc){
			if($item =~ /(.+)\(.+\)/){
				$entries{$1} = 1;
			}else{
				next;
			}				
		}
		foreach my $item (sort keys %entries){
			push(@entries, $item);	
		}		
	 	$out_ncenc = join(',',@entries);
	}else{
		$out_ncenc = '';
	}

	#hot-----------------------------------------------------
	if($info_hot){
		$out_hot = 'Y';	
	}else{
		$out_hot = '';
	}
	#GENE=FAM200B(Distal)[pearson(H3K27ac):0.918602,pearson(H3K4me1):0.912105],FBXL5(Distal)[pearson(H3K27ac):0.875722,pearson(H3K4me1):0.873539]
	#gene-----------------------------------------------------
	if($info_gene){
		my %entries; my @entries;
		my ($header, $data) = split('=',$info_gene);
		my @info_gene = split('],', $data);
		
		foreach my $item (@info_gene){
			if($item =~ /(.+)\((.+)\)\[.*/){
				my $gene = $1 . '(' . $2 . ')';
				$entries{$gene} = 1;
			}elsif($item =~ /(.+)\((.+)\)/){
				my $gene = $1 . '(' . $2 . ')';
				$entries{$gene} = 1;
			}else{
				next;
			}				
		}
		foreach my $item (sort keys %entries){
			push(@entries, $item);	
		}		
	 	$out_gene = join(',',@entries);		
	}else{
		$out_gene = '';
	}
	
	#recur-----------------------------------------------------
	if($info_recur){
		$out_recur = 'Y';
	}else{
		$out_recur = '';
	}
	my $out_string = $out_gerp  . "\t" . 
					 $out_hub   . "\t" .
					 $out_ncenc . "\t" .
					 $out_hot   . "\t" .
					 $out_gene  . "\t" .
					 $out_recur;

	return $out_string;
}