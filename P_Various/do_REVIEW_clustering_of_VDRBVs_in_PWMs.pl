#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
#use Algorithm::Combinatorics qw(combinations);

#example use of combinations
#chdir $dir;
#my @files = <*.*>;
#my $combinations = combinations(\@files,2);
#while (my $pair = $combinations->next) {

#1st QUESTION - how many VDR-BVs in each RXRA:VDR motif?
#I want to plot a histogram depicting how many intersections with a VDR-BV does every RXRA:VDR motif instance have.
#So, I know only a proportion of VDR-BV intersect strong motifs at all. 
#Are there any such motifs hit by MORE than 1 VDR-BVs? If so what's their score? Where are they?

#2nd QUESTION:how many VDR-BV hit 1)no PWM 2)one PWM ... 3)n DISTINCT PWMs?

my $BEDTOOLS = `which bedtools`; chomp $BEDTOOLS;
my $RSCRIPT = `which RRscript`; chomp $RSCRIPT;

my $IN_VDRBV_BED = '/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Output_noDBRECUR.bed';
my $IN_VDRrBV_BED = '/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Output_VDR_rBVs_hg19.bed';
my $PWM_FILE = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Processed_PFMs_jaspar_FUNSEQ_INPUT.txt";
my $INPUT_PSCANCHIP_RIS_DIR;
my $identifier;
my $motif_name;
my $INPUT_VAR_BINARY;
my $MIN_SCORE;
my $input_variants;

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
}elsif($INPUT_VAR_BINARY eq 'rBV'){
	$input_variants = $IN_VDRrBV_BED;
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

####################
#1 get all PWM intervals from the directory, check length against pwms
####################
my %pwm_intervals; #pwm_intervals{PWM_i}{interval_j} = 1
chdir $INPUT_PSCANCHIP_RIS_DIR;
my @files = <Pscanchip*.ris>;
foreach my $FILE (@files){
	print $FILE, "\n";
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
	print STDERR "The length of the motif: $full_motif_id according to the JASPAR Pwm is $MOTIF_LENGTH\n";	
	
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


#make a small histogram for each pwm showing the number of VDR-BVs per motif interval
foreach my $this_pwm (sort keys %pwm_intervals){
	print $this_pwm, "\n";
	my $temp_pwm_bed   = $INPUT_PSCANCHIP_RIS_DIR . '/TEMP_' . $this_pwm . '.bed';
	my $temp_Rcode     = $INPUT_PSCANCHIP_RIS_DIR . '/TEMP_' . $this_pwm . '.R';
	my $temp_Rdata     = $INPUT_PSCANCHIP_RIS_DIR . '/TEMP_' . $this_pwm . '.Rdata';	
	my $out_Rplot      = $INPUT_PSCANCHIP_RIS_DIR . '/plot_vdrbvhist_' . $this_pwm . '.png';
	open (my $outstream,  q{>}, $temp_pwm_bed) or die("Unable to open $temp_pwm_bed : $!");
	foreach my $this_pwm_interval (sort keys %{$pwm_intervals{$this_pwm}}){
		print $outstream $this_pwm_interval, "\n";
	}
	close $outstream;
	system "$BEDTOOLS intersect -c -a $temp_pwm_bed -b $input_variants | cut -f 4 | sort > $temp_Rdata";
	#write R code 
	open ($outstream,  q{>}, $temp_Rcode) or die("Unable to open $temp_Rcode : $!");
	print $outstream "data <- read.table(\"$temp_Rdata\",sep=\"\\t\")" . "\n";
	print $outstream "sub_data <- subset(data, V1 > 0)" . "\n";
	print $outstream "png(file=\"$out_Rplot\")" . "\n";
	print $outstream "hist(sub_data\$V1, xlab=\"#INPUT_VAR_BINARY per motif instance\", ylab=\"Motif Count\", col=\"blue\", main=\"$INPUT_VAR_BINARY cardinality in $this_pwm instances\")" . "\n";
	print $outstream "dev.off()" . "\n";
	close $outstream;
	system "$RSCRIPT $temp_Rcode";
	
	unlink $temp_pwm_bed;
	unlink $temp_Rdata;
	unlink $temp_Rcode;
	#exit;
}
