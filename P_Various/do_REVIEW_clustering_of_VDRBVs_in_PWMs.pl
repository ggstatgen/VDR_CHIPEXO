#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;

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
my $PLOT_EXT = 'png';

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

	#process vcf lines and save in structure (for point 4)
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


#####
##4 go through VDR-BVs and print tsv. with three additional fields IF there is at least 1 intersection. Print: 
#a. number of intersections with unique enriched PWMs
#b. list of PWM ids
#c. whether RXRA::VDR is one of them
#####
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#I want chrom, pos, id, ref, alt, INFO(gene) + new ones (n.pwms hit?which pwms?is vdr one of them?)
my $temp_out_tsv = $INPUT_PSCANCHIP_RIS_DIR . '/' .  $INPUT_VAR_BINARY . '_hittingPWMs_minPWMscore_' . $MIN_SCORE . '_temp.data';
my $out_Rdata    = $INPUT_PSCANCHIP_RIS_DIR . '/' .  $INPUT_VAR_BINARY . '_hittingPWMs_minPWMscore_' . $MIN_SCORE . '.Rdata';
my $out_Rcode    = $INPUT_PSCANCHIP_RIS_DIR . '/' .  $INPUT_VAR_BINARY . '_hittingPWMs_minPWMscore_' . $MIN_SCORE . '.R';
my $out_Rplot    = $INPUT_PSCANCHIP_RIS_DIR . '/' .  $INPUT_VAR_BINARY . '_hittingPWMs_minPWMscore_' . $MIN_SCORE . '.png';
my $out_tsv_nH   = $INPUT_PSCANCHIP_RIS_DIR . '/' .  $INPUT_VAR_BINARY . '_hittingPWMs_minPWMscore_' . $MIN_SCORE . 'nH.tsv';
my $out_tsv      = $INPUT_PSCANCHIP_RIS_DIR . '/' .  $INPUT_VAR_BINARY . '_hittingPWMs_minPWMscore_' . $MIN_SCORE . '.tsv';
my $HEADER = "#CHROM\tPOS\tID\tREF\tALT\tGENE_INFO\tN_PWM_HIT\tPWM_IDs\tVDR_DR3_PRESENT\tCTCF_PRESENT\tZnf423_PRESENT\n";;
open (my $outstream,  q{>}, $temp_out_tsv) or die("Unable to open $temp_out_tsv : $!");
open ($instream,  q{<}, $input_variants_vcf) or die("Unable to open $input_variants_vcf : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/);
	next if($_ eq '');
	my @line = split("\t",$_);
	my $coord = $line[0] . '-' . $line[1];
	my @pwmstring; my $pwmlist = ''; my $counter = 0; 
	my $VDR_DR3_FOUND = '';
	my $CTFC_FOUND = '';
	my $ZNF423_FOUND = '';
	
	#now search the catalogue of intersections with this line
	if($VDRBV_intersections{$coord}){
		foreach my $pwm (keys %{ $VDRBV_intersections{$coord} }){
			if($pwm =~ /MA0074.1/){ $VDR_DR3_FOUND = 'Y'; } 
			if($pwm =~ /MA0116.1/){ $ZNF423_FOUND = 'Y'; }
			if($pwm =~ /MA0139.1/){ $CTFC_FOUND = 'Y'; }					
			push(@pwmstring,$pwm);
			$counter++; 
		}
		my @sorted_pwmstring = sort @pwmstring;
		$pwmlist = join(',',@sorted_pwmstring);																
	}
	
	my $gene_info = get_gene_info($line[7]);
	my $line =  $line[0]   . "\t" . 
				$line[1]   . "\t" . 
				$line[2]   . "\t" .
				$line[3]   . "\t" .
				$line[4]   . "\t" .
				$gene_info . "\t" .
				$counter   . "\t" .
				$pwmlist   . "\t" . 
				$VDR_DR3_FOUND . "\t" .
				$CTFC_FOUND    . "\t" .
				$ZNF423_FOUND;	
	print $outstream $line, "\n";
}
close $instream;
close $outstream;

#sort and get unique out.
system "cat $temp_out_tsv | sort -k1,1V -k2,2n | uniq > $out_tsv_nH";
unlink $temp_out_tsv;

#Write R bar plots
system "cat $out_tsv_nH | cut -f 7 | sort > $out_Rdata";
open ($outstream,  q{>}, $out_Rcode) or die("Unable to open $out_Rcode : $!");
print $outstream "data <- read.table(\"$out_Rdata\",sep=\"\\t\")" . "\n";
print $outstream "data_c <- table(data)" . "\n";
print $outstream "$PLOT_EXT(file=\"$out_Rplot\")" . "\n";
print $outstream "bplt <- barplot(data_c, xlab=\"\#PWM models hit by $INPUT_VAR_BINARY\", ylab=\"\#$INPUT_VAR_BINARY\", width=1, main=\"\#$INPUT_VAR_BINARY by \#PWM models hit(t=$MIN_SCORE)\")" . "\n";
print $outstream "text(x=bplt, y=data_c, labels=as.character(data_c), pos=3, cex = 0.8, col = \"red\", xpd=TRUE)" . "\n";
print $outstream "dev.off()" . "\n";
close $outstream;

system "$RSCRIPT $out_Rcode";
system "sed 1i\"$HEADER\" $out_tsv_nH > $out_tsv";

unlink $out_tsv_nH;
unlink $temp_out_tsv;
unlink $out_Rcode;
unlink $out_Rdata;



#subs-------------------------------------------------------------------
#only keep GENE info, if available, from the funseq2 annotation
#SAMPLE=interestingHets_NA06986_EBLfiltered_hg19;GERP=1.01;CDS=No;HUB=EXOSC7:PHOS(0.947)PPI(0.739);GENE=EXOSC7(Intron);NCDS=0.844214221203786
sub get_gene_info{
	my ($this_info_string) = @_;
	my $gene_field = 'NA';
	
	my @fields = split(";", $this_info_string);
	foreach my $item (@fields){
		$gene_field = $item if($item =~ /^GENE/);
	}
	return $gene_field;
}
