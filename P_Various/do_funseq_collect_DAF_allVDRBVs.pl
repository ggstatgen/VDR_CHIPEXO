#!/usr/bin/perl
use strict;
use warnings;
use List::MoreUtils qw(uniq);
use Getopt::Long;
use File::Basename;
use IO::Zlib;

#25/5/2015 Chris wants me to look at the DAF for all VDR-BVs/VDR-hcBVs not in the motif. I will discard those in the VDR-RXR motif

#23/3/2015 Chris wants me to discard all positions that do not have ancestral state. Maybe count them to see how many are there?

#25/2/2015 
#I collected and plotted phastcons scores for each nucleotide at the motif positions.
#This gives me a perspective over vertebrate conservation at the motif

#I now want to see what happens at the population level: I will collect derived allele frequencies for the VDR-BV variants.
#The hypothesis is that LOB and GOB, outside of the VDR:RXR motif, don't show significant difference, unlike GOB and LOB inside the motif

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
#or plot it as "inconclusive"?

#http://www.1000genomes.org/category/frequently-asked-questions/population
#CEU are EUR, YRI are AFR

#input: funseq file, variants tested 

#algorithm
#-save 1kg variants in hash. Data to keep is AFR_AF and EUR_AF

my $infile;
my $INPUT_VARIANTS = '/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/ALLELESEQ_20samples_merged.vcf.gz'; #b37
my $INPUT_ALLELESEQ_RAW = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/_RAW/interestingHets_vdrchipexo.txt"; #b37
my $ENH_ONLY;
my $subset_bed; #if you only want to collect data for positions in an external bed file (eg CLASS I sites)
my $sym_variants; #flag, set it if you're dealing with sym variants

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
	$output  = $directory . $basename . '_DAF_noDR3motif_subsetbed.Rdata';
}elsif($ENH_ONLY){
	$output  = $directory . $basename . '_DAF_noDR3motif_enhonly.Rdata';
}else{
	$output  = $directory . $basename . '_DAF_noDR3motif.Rdata';
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
#1 -build  hash to map chr-pos to DAF-CEU and DAF-YRI
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
	my $ALTF_AFR; my $ALTF_EUR; #ALTERNATE allele frequencies
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
	if($info[0] =~ /^AA=(.*)/){
		my $anc = uc($1);
		if($anc eq $alt){
			$FLAG = 1;
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
	#get allele frequencies========================
	$ALTF_AFR = $1 if($info =~ /AFR_AF=([0-9]+\.[0-9]+)/);
	$ALTF_EUR = $1 if($info =~ /EUR_AF=([0-9]+\.[0-9]+)/);	
	
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
}
close FILE;

###################
#3 get affinity binding phenotype info 
###################
#chr-position => samplename => alleles / read counts info
my %position2sample2readdepth;
open (my $instream,  q{<}, $INPUT_ALLELESEQ_RAW) or die("Unable to open $INPUT_ALLELESEQ_RAW : $!");
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

#for each coordinate, there might be multiple phenotypes.
#you want to attach one phenotype label per key
my %data_VDRBV2DAF;
open ($instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/);
	next if($_ eq '');
	#I'm looking for rows that DO NOT HAVE THE FOLLOWING:
	#NCENC=TFM(VDR|VDR_JASPAR|chr7:5013558-5013573),TFM(VDR|VDR_XXmotif|chr7:5013558-5013573)...VDR_dreme
	#if there is no such line, there is no motifbr for vdr either
	next if($_ =~ /TFM\(VDR\|VDR/);

	my $LABEL_DIR; my $LABEL_RECUR = 'N';	
	my $info_ncenc;	my $info_motifbr;
	my ($chr, $pos, $rs_id, $info) = (split /\t/)[0,1,2,7]; #chr is hg19
	#FILTER: BED
	if($subset_bed){ next unless(defined check_coords_in_bed($chr, $pos)); }
	
	#$chr =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19
	my $genomic_coord = $chr . '-' . $pos;
	my @info = split(";", $info);
	
	#I want to flag rows which are in enhancers or break a motif (any but the VDRs)
	foreach my $item (@info){
		$info_ncenc = $item if($item =~ /^NCENC/);
		$info_motifbr = $item if($item =~ /^MOTIFBR/);	
		$LABEL_RECUR = 'Y' if($item =~ /RECUR/);
	}
	#FILTER: ENHANCER
	if($ENH_ONLY){ next unless($info_ncenc =~ /Enhancer/); }
	
	#if you have motifbr data,you can remove instances where motif pwm change is discordant with phenotype change.
	#if you do not have motifbr data, you cannot do this, and you will have to look at the final results
	#
	#motifbr maybe you can do some selection based on vdr-hcBVs which DO NOT impact a motif (any)
#	my ($ncenc,$ncenc_data) = split("=", $info_ncenc);
#	my @ncenc_data = split(",", $ncenc_data);
	
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
		$LABEL_DIR = 'LOB';
	}elsif( ($gob_samples{$genomic_coord} > 0 ) && ($lob_samples{$genomic_coord} == 0)  ){
		$LABEL_DIR = 'GOB';
	}else{
		print STDERR "ERROR VDR-BV $rs_id $lob_samples{$genomic_coord} lob samples and  $gob_samples{$genomic_coord} gob samples - should not be here.\n";
		next;
	}
						
	if($variants_1kg{$genomic_coord}{DAFEUR}){
		$data_VDRBV2DAF{$genomic_coord}{DAFEUR} = $variants_1kg{$genomic_coord}{DAFEUR};
	}else{
		$data_VDRBV2DAF{$genomic_coord}{DAFEUR} = 'NA';
	}
	if($variants_1kg{$genomic_coord}{DAFAFR}){
		$data_VDRBV2DAF{$genomic_coord}{DAFAFR} = $variants_1kg{$genomic_coord}{DAFAFR};
	}else{
		$data_VDRBV2DAF{$genomic_coord}{DAFAFR} = 'NA';
	}
	#save labels
	$data_VDRBV2DAF{$genomic_coord}{LABEL_DIR} = $LABEL_DIR;
	$data_VDRBV2DAF{$genomic_coord}{LABEL_RECUR} = $LABEL_RECUR;

}
close $instream;

#motifpos -> genomic_coord -> DAFCEU,...DAFASN
#motifpos -> genomic_coord -> LABEL_DIR

my $TYPE;
if($sym_variants){
	$TYPE = 'sym';
}else{
	$TYPE = 'asym';
}

open (my $outstream,  q{>}, $output) or die("Unable to open $output : $!");
print $outstream "TYPE\tDIRECTION\tETHNICITY\tDAF\tIS_VDRhcBV\n";

foreach my $genomic_pos (sort keys %data_VDRBV2DAF ){
	my $direction = $data_VDRBV2DAF{$genomic_pos}{LABEL_DIR};
	my $is_vdrhcbv = $data_VDRBV2DAF{$genomic_pos}{LABEL_RECUR};
        my $eur = $data_VDRBV2DAF{$genomic_pos}{DAFEUR};
        my $afr = $data_VDRBV2DAF{$genomic_pos}{DAFAFR};
        print $outstream $TYPE . "\t" . $direction . "\t" . 'CEU' . "\t" . $eur . "\t" .  $is_vdrhcbv . "\n";
      	print $outstream $TYPE . "\t" . $direction . "\t" . 'YRI' . "\t" . $afr . "\t" .  $is_vdrhcbv . "\n";
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
