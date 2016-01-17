#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#16/1
#Now I incorporate code to do both checks in here:
#check pwm in interval actually hits the vdrbv
#check pwm in interval has score >= optional score threshold

#NOTE 
#
#6/1/2016
#Having shown the results to Chris, he says they are too noisy. 
#So I want to implement a new feature: I want this to be able to output only intervals within a certain pscanchip score

#/11/12/2015
#I FOUND OUT THAT THE INTERVAL IN THE PSCANCHIP ris FILE CAN BE WRONG, DEPENDING ON THE TF considered!!!
##The only way to solve this is to have a file with the ACTUAL pwm sizes from Jaspar, indications are on the google drive doc
#I generated one from the Jaspar and I link it here. You need to slurp it and count the sizes of the PWM and check the interval in the .ris file is correct
#OTHERWISE FUNSEQ2 complains

#quick script to turn the Pscanchip .ris motif files in a format compatible with the funseq annotation "ENCODE.tf.bound.union.bed"
#create this, then append to the existing annotation in funseq, sort, gzip

#AFTER FAILING TO GET RESULTS from funseq, I noticed that the length of the motif returned from pscanchip is ONE NT SHORTER than the length of the motif inferred by funseq from the PFM. So I now OPEN the start coordinate (-1)
#ATTENTION 20/11: opening the coordinate gives me a truncated motif. I need to ADD ONE TO THE END COORDINATE from pscanchip
##ATTENTION 4/2/2015: it seems I need to add + to both start and end of pscanchip to obtain coordinates which give me the exact sequence in the UCSC browser!!
##ATTENTION 4/2/2015: if I do the latter, I will get a 14 bases motif; also if I do the latter in bedtools get fasta I will get a 14bp motif. Therefore the last way is wrong. Reverting.

#https://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1

#ALSO DO NOT add any "_\d+mer" suffix to the motif name because funseq will chop it and fail to match the motif name

#input format (pscanchip ris)
#This file is in the following format:
#CHR     REG_START       REG_END REG_STRAND      ABS_SITE_START  ABS_SITE_END    REL_SITE_START  REL_SITE_END    SITE_STRAND     SCORE   SITE
#chr5    131986387       131986536       0       131986469       131986483       7       20      -0      0.965192        TGAACCCTGTGACCT
#chr5    1297635 1297784 0       1297697 1297711 -13     0       -0      0.941996        TGAACTCCATGAACT
#chr3    189207271       189207420       0       189207374       189207388       28      41      0       0.9018  AGGCCATTGAGTTCA
#chr6    158182596       158182745       0       158182674       158182688       3       16      -0      0.894597        TGAACCGTGTGACCT
#chr7    139586182       139586331       0       139586226       139586240       -31     -18     -0      0.890353        GGAACCTATTAACCT

#output format
#chr start stop motif name . strand, TF name
#eg
#chr1    714055  714076  ATF3_disc2_8mer .       +       ATF3
#chr1    714064  714078  ATF3_known8_8mer        .       +       ATF3
#chr1    714066  714075  ATF3_known1_8mer        .       +       ATF3

my $infile_vcf;
my $infile_ris;
my $identifier;
my $motif_name;
my $PWM_FILE;
my $MIN_SCORE;
GetOptions(
		'i_vcf=s'	=>\$infile_vcf,    
        'i_ris=s'	=>\$infile_ris,
        'pwm=s'		=>\$PWM_FILE,
        'id=s'		=>\$identifier,
        'm=s'		=>\$motif_name,
        's=f'		=>\$MIN_SCORE
);

#$PWM_FILE="/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Processed_PFMs_jaspar_FUNSEQ_INPUT.txt";

my $USAGE = "USAGE: do_funseq_adapt_motiffile.pl -i_vcf=<INFILE_VDRBV> -i_ris=<INFILE_RIS> -pwm=<ENCODE_PWM_FILE> -id=<ID> -m=<MOTIF_NAME> (opt)-s=<MINSCORE>\n" .
		"<INFILE_VDRBV>  .vcf file out of funseq, with the VDR-BVs"  .
		"<INFILE_RIS> .ris file from PscanChIP or .out file from do_pscanchip_out_intersect_vdrbvs.pl\n" .
		"<ENCODE_PWM_FILE> text file with Jaspar PWMs in ENCODE format obtained with RSAT\n" .
		"<ID> string to use for the ID (e.g. Jaspar ID) of the PWM in the output file\n" . 
		"<MOTIF_NAME> string to use for the Motif PWM name in the output file\n"	.
		"optional <MINSCORE> lower threshold on score (eg 0.8) (default:none)\n";

unless($infile_vcf && $infile_ris && $identifier && $motif_name && $PWM_FILE){
	print $USAGE;
	exit -1;
}
print STDERR "THRESHOLDING ON SCORE: $MIN_SCORE\n" if($MIN_SCORE);

#######
#1 get the REAL motif length from the encode representations of the motif
#######
my %motif; my $A = 1; my $C = 2; my $G = 3; my $T = 4;
my $prev_name; my @info; my $temp;
#slurp pwm file
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


#heterodimers are saved by Jaspar as monomer::monomer
#I replaced the :: with a '-' in the input file name because it's not recognised by the SGE submission
#change again here:
$motif_name =~ s/\-/\:\:/;
#get the length of the motif analyzed in this iteration
my $full_motif_id = $motif_name . '_' . $identifier;
my $ref = $motif{$full_motif_id};
my $MOTIF_LENGTH = scalar(@$ref); 
print STDERR "The length of the motif: $full_motif_id according to the JASPAR Pwm is $MOTIF_LENGTH\n";

################
#2 save the vdr chr-pos in a data structure
################
my %variant_coords;
open ($instream,  q{<}, $infile_vcf) or die("Unable to open $infile_vcf : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/);
	next if($_ eq '');
	
	#get coords
	my ($chr, $pos) = (split /\t/)[0,1];
	unless($chr =~ /^chr/){
		$chr = 'chr' . $chr;
	}
	my $coord = $chr . '-' . $pos;
	$variant_coords{$coord} = 1;
}
close $instream;

############
#3 go through all .ris intervals and only save them if BOTH the following are true:
#a) the interval CONTAINS the variant
#b) the interval has a score >= MIN_SCORE if MIN_SCORE is set
############
my $counter_failscore = 0;
my $counter_failinterval = 0;
open ($instream,  q{<}, $infile_ris) or die("Unable to open $infile_ris : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^CHR/);
	
	my ($chr,$motif_start,$motif_end,$motif_strand,$score,$site) = (split /\t/)[0,4,5,8,9,10];
	#next if (!$chr);
	if(  $MIN_SCORE && ($score < $MIN_SCORE) ){
		$counter_failscore += 1;
		next;
	}
	my $pscanchip_interval_length = ($motif_end - $motif_start);

	if($pscanchip_interval_length <  $MOTIF_LENGTH){
		$motif_end += 1;
	}elsif($pscanchip_interval_length == $MOTIF_LENGTH){	
		;
	}else{
		print STDERR "ERROR: the pscanchip motif length is LARGER than the Jaspar length. Verify.\n";
		exit -1;	
	}
	
	if(interval_contains_vdrbvs($chr, $motif_start, $motif_end, %variant_coords)){
		print STDOUT $chr . "\t" . $motif_start . "\t" . $motif_end . "\t" .  $full_motif_id . "\t" . '.', "\t" . $motif_strand . "\t", $motif_name .  "\n";
	}else{
		$counter_failinterval += 1;
		next;	
	}
}
close $instream;

print STDERR "Number of .ris intervals discarded because of score: $counter_failscore\n";
print STDERR "Number of .ris intervals discarded because no VDR-BV is in the PWM motif: $counter_failinterval\n";
print STDERR 'FINISHED.';

#############
# should return TRUE if the ris interval contains at least a VDR-BV; under otherwise
#############
sub interval_contains_vdrbvs{
	my ($thischr, $int_start, $int_end, %vdrbv_hash) = @_;
	
	foreach my $item (keys %vdrbv_hash){
		my ($vdrbv_chr, $vdrbv_pos) = split('-', $item);
		
		if( ($thischr eq $vdrbv_chr) && ($vdrbv_pos >= $int_start) && ($vdrbv_pos <= $int_end ) ){
			return 1;
		}
	}
	return undef;
}
