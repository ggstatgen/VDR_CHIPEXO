#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
#use Algorithm::Combinatorics qw(combinations);

#example use of combinations
#chdir $dir;
#my @files = <*.*>;
#my $combinations = combinations(\@files,2);
#while (my $pair = $combinations->next) {

#1st QUESTION - how many VDR-BVs in each RXRA:VDR motif?
#I want to plot a histogram depicting how many intersections with a VDR-BV does every PWM motif instance have.
#So, I know only a proportion of VDR-BV intersect strong motifs at all. 
#Are there any such motifs hit by MORE than 1 VDR-BVs? If so what's their score? Where are they?

#2nd QUESTION:how many VDR-BV hit 1)no PWM 2)one PWM ... 3)n DISTINCT PWMs? (like the one before but with all pwms, or a selection?)
#output: annotated tsv, saying for each vdrbv, the list of pwms which are hit and the score.

#for now, save in a tsv. the initial .vcf of variants ONLY when they hit > 1 PWMs

#3rd QUESTION - intorno di vdrbv
#if i pick a vdrbv hitting ANY motif, are there any other good motif PWMs for other tfs in the surroundings??



my $BEDTOOLS = `which bedtools`; chomp $BEDTOOLS;
my $RSCRIPT = `which RRscript`; chomp $RSCRIPT;

my $IN_VDRBV_BED = '/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Output_noDBRECUR.bed';
my $IN_VDRrBV_BED = '/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Output_VDR_rBVs_hg19.bed';
my $IN_VDRBV_VCF = '/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/Output_noDBRECUR.vcf';
my $IN_VDRrBV_VCF = '/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/SUPPL_DATA_Output_noDBRECUR_REP_hg19.vcf';
my $PWM_FILE = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Processed_PFMs_jaspar_FUNSEQ_INPUT.txt";
my $INPUT_PSCANCHIP_RIS_DIR;
my $identifier;
my $motif_name;
my $INPUT_VAR_BINARY;
my $MIN_SCORE;
my $input_variants;
my $input_variants_vcf;

GetOptions(
        'data=s'   		=>\$INPUT_PSCANCHIP_RIS_DIR,        
	'variants=s'		=>\$INPUT_VAR_BINARY,
        'pwmscore=f'		=>\$MIN_SCORE            
);
my $USAGE = "\nUSAGE: $0 -data=<INFILE_PSCANCHIP_DIR> -variants=<BV|rBV> (opt) -pwmscore=<MINSCORE>\n" .
			"<INFILE_PSCANCHIP> directory of .ris files from PscanChip\n" .
			"<BV|rBV> if BV, all VDRBV will be used; if rBV, only VDR-rBV will be used\n" .
			"<MINSCORE> lower threshold on score (eg 0.8) (default:none)\n";
			
unless($INPUT_PSCANCHIP_RIS_DIR && $INPUT_VAR_BINARY){
	print $USAGE;
	exit -1;
}
#binary arguments
if($INPUT_VAR_BINARY eq 'BV'){
	$input_variants = $IN_VDRBV_BED;
	$input_variants_vcf = $IN_VDRBV_VCF;
}elsif($INPUT_VAR_BINARY eq 'rBV'){
	$input_variants = $IN_VDRrBV_BED;
	$input_variants_vcf = $IN_VDRrBV_VCF;
}else{
	print STDERR "ERROR: field -v not recognised: $INPUT_VAR_BINARY. Aborting.\n";
	exit -1;
}
print STDERR "THRESHOLDING ON SCORE: $MIN_SCORE\n" if($MIN_SCORE);
$INPUT_VAR_BINARY = 'VDR-' . $INPUT_VAR_BINARY;

####################
#0 slurp motif lengths
####################
#slurp pwm file
my %motif; my $A = 1; my $C = 2; my $G = 3; my $T = 4;
my $prev_name; my @info; my $temp;	
open (my $instream,  q{<}, $PWM_FILE) or die("Unable to open $PWM_FILE : $!");
while(<$instream>){
	chomp $_;
	if(/^>/){
		$prev_name = (split/>|\s+/,$_)[1];
	}else{
		@info = split/\s+/,$_;
		if(not exists $motif{$prev_name}){
			$motif{$prev_name}->[0] = {(A=>$info[$A], T=>$info[$T], C=>$info[$C], G=>$info[$G])};
		}else{
			$temp = $motif{$prev_name};
			$motif{$prev_name}->[scalar(@$temp)] = {(A=>$info[$A], T=>$info[$T], C=>$info[$C], G=>$info[$G])};
		}
	}
}
close $instream;

#####
##1 slurp VDR-BV vcf indexed by coordinate (needed for question 2)
#####
#my %VDRBV_library;
#open ($instream,  q{<}, $input_variants_vcf) or die("Unable to open $input_variants_vcf : $!");
#while(<$instream>){
#	chomp;
#	next if($_ =~ /^\#/);
#	next if($_ eq '');
#	my @line = split("\t",$_);
#	my $chr = shift @line;
#	my $pos = shift @line;
#	my $line = join("\t", @line);
#	my $coord = $chr . '-' . $pos;
#	$VDRBV_library{$coord}{'VCFDATA'} = $line;
#}


####################
#2 get all PWM intervals from the directory, check length against pwms
####################
my %pwm_intervals; #pwm_intervals{PWM_i}{interval_j} = 1
chdir $INPUT_PSCANCHIP_RIS_DIR;
my @files = <Pscanchip*.ris>;
foreach my $FILE (@files){
	#get motif name and id from the filename
	#eg Pscanchip_hg19_bkgGM12865_Jaspar_VDRBVs_RXRA-VDR_MA0074.1_sites
	if($FILE =~ /BVs_(.*)_(.*)_sites/){
		$motif_name = $1;
		$identifier = $2;
	}else{
		print STDERR "ERROR: Unable to recognise the input motif file name: $FILE. Skipping..\n";
		next;
	}
	#heterodimers are saved by Jaspar as monomer::monomer
	#I replaced the :: with a '-' in the input file name because it's not recognised by the SGE submission
	#change again here:
	$motif_name =~ s/\-/\:\:/;
	#get the length of the motif analyzed in this iteration
	my $full_motif_id = $motif_name . '_' . $identifier;
	my $ref = $motif{$full_motif_id};
	my $MOTIF_LENGTH = scalar(@$ref); 
	
	open ($instream,  q{<}, $FILE) or die("Unable to open $FILE : $!");
	while(<$instream>){
		chomp;
		next if($_ eq '');
		next if($_ =~ /^CHR/);
		my ($chr,$motif_start,$motif_end,$motif_strand,$score,$site) = (split /\t/)[0,4,5,8,9,10];
		next if(  $MIN_SCORE && ($score < $MIN_SCORE) );
		
		my $pscanchip_interval_length = ($motif_end - $motif_start);		
		if($pscanchip_interval_length <  $MOTIF_LENGTH){
			$motif_end += 1;
		}elsif($pscanchip_interval_length == $MOTIF_LENGTH){	
			;
		}else{
			print STDERR "ERROR: the pscanchip motif length is LARGER than the Jaspar length. Verify.\n";
			exit -1;	
		}
		my $this_motif_interval = $chr . "\t" . $motif_start . "\t" . $motif_end;
		$pwm_intervals{$full_motif_id}{$this_motif_interval} = 1;		
	}
	close $instream;
}


########
#3 for each PWM model, save one png R file showing a histogram with the number of VDR-BV hitting each model.
########
unless($MIN_SCORE) { $MIN_SCORE = 'NA'; }
#make a small histogram for each pwm showing the number of VDR-BVs per motif interval
my %VDRBV_intersections;
foreach my $this_pwm (sort keys %pwm_intervals){
	print $this_pwm, "\n";
	my $temp_pwm_bed   = $INPUT_PSCANCHIP_RIS_DIR . '/TEMP_' . $this_pwm . '.bed';
	my $temp_pwm_bed_o = $INPUT_PSCANCHIP_RIS_DIR . '/TEMP_' . $this_pwm . '_out.bed';
	my $temp_Rcode     = $INPUT_PSCANCHIP_RIS_DIR . '/TEMP_' . $this_pwm . '.R';
	my $temp_Rdata     = $INPUT_PSCANCHIP_RIS_DIR . '/TEMP_' . $this_pwm . '.Rdata';	
	my $out_Rplot      = $INPUT_PSCANCHIP_RIS_DIR . '/plot_vdrbvhist_' . $this_pwm . '_' .  $INPUT_VAR_BINARY . '_minPWMscore_' . $MIN_SCORE . '.png';
	open (my $outstream,  q{>}, $temp_pwm_bed) or die("Unable to open $temp_pwm_bed : $!");
	foreach my $this_pwm_interval (sort keys %{$pwm_intervals{$this_pwm}}){
		print $outstream $this_pwm_interval, "\n";
	}
	close $outstream;
	system "$BEDTOOLS intersect -a $input_variants_vcf -b  $temp_pwm_bed > $temp_pwm_bed_o";
	system "$BEDTOOLS intersect -c -a $temp_pwm_bed -b $input_variants | cut -f 4 | sort > $temp_Rdata";

	#process vcf lines and save in structure
	open (my $instream,  q{<}, $temp_pwm_bed_o) or die("Unable to open $temp_pwm_bed_o : $!");
	while(<$instream>){
		chomp $_;
		my @line = split("\t",$_);
		my $chr = shift @line;
		my $pos = shift @line;
		my $coord = $chr . '-' . $pos;
		$VDRBV_intersections{$coord}{$this_pwm} = 1;
	}
	close $instream;


	#write R code 
	open ($outstream,  q{>}, $temp_Rcode) or die("Unable to open $temp_Rcode : $!");
	print $outstream "data <- read.table(\"$temp_Rdata\",sep=\"\\t\")" . "\n";
#	print $outstream "sub_data <- subset(data, V1 > 0)" . "\n";
	print $outstream "data_c <- table(data)" . "\n";
	print $outstream "png(file=\"$out_Rplot\")" . "\n";
	print $outstream "bplt <- barplot(data_c, xlab=\"\#$INPUT_VAR_BINARY per motif instance\", ylab=\"Motif Count\", width=1, main=\"$INPUT_VAR_BINARY cardinality ($this_pwm; thrs=$MIN_SCORE)\")" . "\n";
	print $outstream "text(x=bplt, y=data_c, labels=as.character(data_c), pos=3, cex = 0.8, col = \"red\", xpd=TRUE)" . "\n";
	print $outstream "dev.off()" . "\n";
	close $outstream;
	system "$RSCRIPT $temp_Rcode";
	
	unlink $temp_pwm_bed;
	unlink $temp_Rdata;
	unlink $temp_Rcode;
	unlink $temp_pwm_bed_o;
}

#plot .tsv of vdr-bv positions followed by a string indicating all the UNIQUE PWMs they hit
my $out_tsv = $INPUT_PSCANCHIP_RIS_DIR . '/' .  $INPUT_VAR_BINARY . '_hittingPWMs_minPWMscore_' . $MIN_SCORE . '.tsv';
open (my $outstream,  q{>}, $out_tsv) or die("Unable to open $out_tsv : $!");
print $outstream "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tPWM_HIT\n";
foreach my $vdr_coord (sort  { $a <=> $b } keys %VDRBV_intersections){
	my @pwmstring; my $pwmlist;
	foreach my $pwm (keys %{ $VDRBV_intersections{$vdr_coord} }){ push(@pwmstring,$pwm); }
	$pwmlist = join(',',@pwmstring);
	my ($chr, $pos) = split('-',$vdr_coord);
	print $outstream $chr . "\t" . $pos . "\t" . $pwmlist, "\n";
}
close $outstream;

#print Dumper(\%VDRBV_library);
#exit;
####
#3 print out a tsv file of those vdr-bv falling in at least 1 DISTINCT PWM model
####
#my $out_tsv      = $INPUT_PSCANCHIP_RIS_DIR . '/' .  $INPUT_VAR_BINARY . '_hittingPWMs_minPWMscore_' . $MIN_SCORE . '.tsv';
#open (my $outstream,  q{>}, $out_tsv) or die("Unable to open $out_tsv : $!");
#print $outstream "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tPWM_HIT\n";
#foreach my $item (keys %VDRBV_library){
#	my @pwmlist; #list of comma separated PWMs that are hit
#	my $pwmstring;
#	foreach my $pwm_hit (sort keys %{$VDRBV_library{$item}{'PWM'}}){ push (@pwmlist,$pwm_hit);  }
#	$pwmstring = join(",", @pwmlist);
#	next if(!$pwmstring || ($pwmstring = '') );
#	
#	my ($chr, $pos) = split('-', $item);
#	my $string = $chr . "\t" . $pos . "\t" . $VDRBV_library{$item}{'VCFDATA'} . "\t" . $pwmstring;
#	print $outstream $string, "\n";
#}
#close $outstream;

#print output
#my $out_table = $INPUT_PSCANCHIP_RIS_DIR . '/TABLE_vdrbv_cardinality.tsv';
#open (my $outstream,  q{>}, $out_table) or die("Unable to open $out_table : $!");
#foreach my $item (%tsv_line){ print $outstream $item, "\n"};
#close $outstream;
