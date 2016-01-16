#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;


#When you have built a table using do_funseq_collect_LOB_GOB_DAF_MOTIFBR.pl, you might want to create a latex table using specific filters on:
#-DAF CEU (if the sample where the VDR-BV is is CEU)
#-DAF YRI (if the sample where the VDR-BV is is YRI)
#motif score

#You could for example want to have a table in the paper with all the VDR-BV with very low DAF in either CEU or YRI and hitting a very conserved motif

my $infile;
my $MOTIF_THRESHOLD = 0.8;
my $DAF_THRESHOLD = 0.1; #if lob, this is the threshold, otherwise if gob threshold = 1 - this
GetOptions(
        'i=s'      =>\$infile
);

#CHR     POS     ID      REF     ALT     SAMPLE  MOTIF_BREAK     MOTIF_STRAND    MOTIF_SCORE     MOTIF_SEQ(+)    MOTIF_HIT_POS   RD_REF_ALT      PH_DIR  PH_FC   DAF_EUR DAF_YRI 1KG_ANNOTATION  FUNSEQ_GERP     FUNSEQ_PROTEIN_NET_HUB  FUNSEQ_NON_CODING_ANNOTATION    FUNSEQ_HOT_REG  FUNSEQ_GENE     FUNSEQ_RECUR

#chr15   93461371        rs1406714       G       C       NA19249(YRI)    N       +       0.671401        GGGAAACCGCGATTA 8       9-1     LOB     5       0.41    0.54    AA=G;AF=0.64;AFR_AF=0.46;AMR_AF=0.65;ASN_AF=0.84;AVGPOST=0.9968;ERATE=0.0005;EUR_AF=0.59;LDAF=0.6367;RSQ=0.9951;SNPSOURCE=LOWCOV;THETA=0.0002;VT=SNP;AN=40;AC=23        -5.98   CHD2    Enhancer,TFM,TFP        Y       CHD2(Intron) 

if(!$infile){
	print "USAGE: do_funseq_collect_LOB_GOB_DAF_MOTIFBR_build_table.pl -i=<INFILE>\n";
    print "<INFILE> output table of do_funseq_collect_LOB_GOB_DAF_MOTIFBR.pl\n";
    print "EDIT FILE TO CHANGE THRESHOLDS. MOTIF: $MOTIF_THRESHOLD; DAF: $DAF_THRESHOLD\n";
    print "discarding non motifbr lob and skipping variants in pos 7,8,9. Edit script to change\n";
    exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	if($_ =~ /^CHR/){ print $_, "\n"; next; }
	my $DAF_THRS;
	
	my @fields = split("\t", $_);
	
	#skip variants hitting the dr3 sequence
	next if($fields[10] eq '7'  || $fields[10] eq '8' || $fields[10] eq '9' );
	
	
	#check if gob or lob. If lob, check it's motifbr. if gob, change threshold
	if($fields[12] =~ /LOB/){
		next unless ($fields[6] =~ /Y/); #must be a motif break
		$DAF_THRS = $DAF_THRESHOLD;
	}elsif($fields[12] =~ /GOB/){
		$DAF_THRS = (1 - $DAF_THRESHOLD);
	}else{
		print STDERR "ERROR: neither lob nor gob?? Shouldn't be here. Aborting..\n";
		exit -1;
	}
	
	#get sample ethnicity and check DAF
	if($fields[5] =~ /YRI/){
		next if(!$fields[15]);
		if($fields[12] =~ /GOB/){
			next unless($fields[15] >= $DAF_THRS);
		}else{
			next unless($fields[15] <= $DAF_THRS);			
		}
	}elsif($fields[5] =~ /CEU/){
		next if(!$fields[14]);
		if($fields[12] =~ /GOB/){
			next unless($fields[14] >= $DAF_THRS);
		}else{
			next unless($fields[14] <= $DAF_THRS);			
		}	
	}else{
		print STDERR "ERROR: no ethnicity available for sample: $fields[5]. Aborting..\n";
		exit -1;
	}
	#check motif score
	next if(!$fields[8]);
	next unless ($fields[8] >= $MOTIF_THRESHOLD);
	
	print $_, "\n";
}
close $instream;