#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use IO::Zlib;


#for the 7 pics intersection, I want the ancestral, the daf and the directionality of binding
#I will further process the latex table obtained with do_funseq_BROAD_PICS_get_latex_table.pl

#ancestral info:
my $INPUT_VARIANTS = '/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/ALLELESEQ_20samples_merged.vcf.gz'; #b37
#read depth info:
my $INPUT_ALLELESEQ_RAW = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/_RAW/interestingHets_vdrchipexo.txt"; #b37
my $LOB = 'lob'; #ancestral -> alternative reduces binding
my $GOB = 'gob'; #ancestral -> alternative increase binding


my $infile;
GetOptions(
        'i=s'      =>\$infile
);
#PICS info:
#recur
#$infile = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/vdrbv_pics/Output_recur_INTERSECT_BROAD_PICS.latex.txt"; #hg19
#all
$infile = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/vdrbv_pics/Output_sorted_INTERSECT_BROAD_PICS_2.tsv"; #hg19

if(!$infile){
	print "USAGE: do_funseq_collect_PICS_DAF_LOBGOB.pl -i=<INFILE>\n";
    print "<INFILE> latex table with intersection of funseq vcf and PICS list\n";
    exit 1;
}
my ($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
my $output  = $directory . $basename  . '_DAF_LOBGOB.latex.txt';

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
	my $ALTF_EUR; my $ALTF_AFR; #ALTERNATE allele frequencies
	my $FLAG; # set to one if the ancestral is the alternate	
	
	my @fields = split("\t", $_);
	#skip indels
	next if(length($fields[3]) > 1);
	next if(length($fields[4]) > 1);
	
	my $key = $fields[0] . '-' . $fields[1];
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

#chr-position => samplename => alleles / read counts info
my %position2sample2readdepth;
open (my $instream,  q{<}, $INPUT_ALLELESEQ_RAW) or die("Unable to open $INPUT_ALLELESEQ_RAW : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^sample/); #header
	next if($_ eq '');
	my ($sample_id, $chr,$snppos, $ref, $m_allele, $p_allele, $cA, $cC, $cG, $cT, $winner, $symcls) = (split /\t/)[0,1,2,3,8,9,10,11,12,13,14,15];
	next if(!$chr);
	next if(!$snppos);
	next unless ($symcls =~ /Asym/);
	
	my $coordinate_id = $chr . '-' . $snppos;
	
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

#open latex table file with VDR-hcBVs which are pics

#CHR
#POS
#RS_ID
#FUNSEQ_ANN
#FUNSEQ_SCORE
#PICS_SCORE
#PICS_TRAIT
#PICS_ANN


#chr1
#180950869
#rs2331903
#GERP=-5.2;HUB=STX6;NCENC=DHS,Enhancer,TFP;HOT=Y;GENE=STX6;RECUR=Y
#2.260
#0.048
#Progressive_supranuclear_palsy
#RISK_ALLELE=T;ANNOTATION=none;NEAREST_GENE=GPR65;EQTL=none;EQTL_DIR=NA;TOP_ENH=CD25-_IL17-_Th_stim_MACS;GM12878

my %data_entries;
open ($instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\CHR/);
	next if($_ eq '');
	my %lob_samples; 
	my %gob_samples; 
	my $BINDING_DIR = '-';my $DAF_CEU = '-'; my $DAF_AFR = '-';
	
	
	my ($chr, $pos, $rs_id, $funseq_ann, $funseq_score, $pics_score, $pics_trait, $pics_ann) = (split /\t/)[0,1,2,3,4,5,6,7];
	my $chr_b37 = $chr; 
	$chr_b37 =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19	
	my $genomic_coord = $chr_b37 . '-' . $pos;
	if(!$position2sample2readdepth{$genomic_coord}){
		print STDERR "$genomic_coord - $rs_id: no ancestral/derived info available for this VDR-BV...Saving as is\n";
		my $line = $_ . "\t" . $BINDING_DIR . "\t" . $DAF_CEU . "\t" . $DAF_AFR;
		$data_entries{$line} = 1;		
		next;
	}
	#############
	#get LOB GOB
	#############
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
	#############
	#decide whether to keep or discard
	#############
	if(  ($lob_samples{$genomic_coord} > 0 ) && ($gob_samples{$genomic_coord} > 0)  ){
		print STDERR "Warning VDR-BV $rs_id - at this position, $lob_samples{$genomic_coord} sample are LOB and $gob_samples{$genomic_coord} are GOB.\n";
		if($lob_samples{$genomic_coord}/$gob_samples{$genomic_coord}  >= 3){
			print STDERR "Choosing LOB.\n";
			$BINDING_DIR = $LOB;
		}elsif($gob_samples{$genomic_coord}/$lob_samples{$genomic_coord}  >= 3){
			print STDERR "Choosing GOB.\n";
			$BINDING_DIR = $GOB;
		}else{
			#print STDERR "Skipping\n";	
			$BINDING_DIR = '-';					
		}
	}elsif( ($lob_samples{$genomic_coord} > 0 ) && ($gob_samples{$genomic_coord} == 0) ){
		$BINDING_DIR = $LOB;
	}elsif( ($gob_samples{$genomic_coord} > 0 ) && ($lob_samples{$genomic_coord} == 0)  ){
		$BINDING_DIR = $GOB;
	}else{
		print STDERR "ERROR VDR-BV $rs_id $lob_samples{$genomic_coord} lob samples and  $gob_samples{$genomic_coord} gob samples - should not be here.\n";
		next;
	}
	############
	#get daf
	############
	$DAF_CEU = $variants_1kg{$genomic_coord}{DAFEUR} if($variants_1kg{$genomic_coord}{DAFEUR});
	$DAF_AFR = $variants_1kg{$genomic_coord}{DAFAFR} if($variants_1kg{$genomic_coord}{DAFAFR});
	
	############
	#build new output
	############
	my $line = $_ . "\t" . $BINDING_DIR . "\t" . $DAF_CEU . "\t" . $DAF_AFR;
	$data_entries{$line} = 1;
}
close $instream;


open (my $outstream,  q{>}, $output) or die("Unable to open $output : $!");
print $outstream "CHR\tPOS\tRS_ID\tFUNSEQ_ANN\tFUNSEQ_SCORE\tPICS_SCORE\tPICS_TRAIT\tPICS_ANN\tBINDING_DIR\tDAF_CEU\tDAF_YRI\n";
foreach my $item (sort keys %data_entries){
	print $outstream $item, "\n";
}
close $outstream;



sub process_read_data{
	my ($read_data_string, $lcl_sample, $id) = @_;
	
	my ($ancestral_ref,$ancestral_alt,$cA,$cC,$cG,$cT) = split(",", $read_data_string);
	my $anc_count; my $der_count;
	my $der_base = $ancestral_alt;
	my $anc_base = $ancestral_ref;
	
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