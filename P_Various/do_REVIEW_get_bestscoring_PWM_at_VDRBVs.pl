#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#6/4/2016
#Following the latest batch of comments from Chris and Wilfried, I want to :
#1-associate each VDR-BV hitting a VDR-RXR motif with strong score, and move aside
#2-for all remaining VDR-BVs, what is the strongest motif (by score) being hit? associate to that

#in the end you want a plot mapping every VDR-BV to the BEST motif they hit.
#All the bottom code is useless atm. Need to rewrite.

#Inputs: all the enriched VDR-BV (not VDR-rBV).ris files with a pwm interval for each of the 43.000
#main loop: for each VDR-BV:
#does the VDR-BV intersect a RXR::VDR with pscanchip > 0.8?
#yes: store in hash with labels
#no: check all ris files. For all the pwms for which an intersection is found, check score (KEEP THRESHOLD AT 0.8?)
#store best scoring  in hash with labels 


my $BEDTOOLS = `which bedtools`; chomp $BEDTOOLS;
my $RSCRIPT = `which RRscript`; chomp $RSCRIPT;

my $CHROMSIZES = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom_simple.sizes';
my $IN_VDRBV = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/Output_noDBRECUR.vcf"; 
#my $IN_VDRBV_BED = '/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Output_noDBRECUR.bed';
#my $IN_VDRrBV_BED = '/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Output_VDR_rBVs_hg19.bed';
my $PWM_FILE = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/Processed_PFMs_jaspar_FUNSEQ_INPUT.txt";
my $RXR_VDR_RIS = "Pscanchip_hg19_bkgGM12865_Jaspar_VDRBVs_RXRA-VDR_MA0074.1_sites.ris";

my $INPUT_RIS_DIR;
my $MIN_SCORE;

#forse rimuovi
my $temp_vars;
my $input_variants;

my $identifier;
my $motif_name;
my $motif_string;

GetOptions(
        'm=s'		=>\$INPUT_RIS_DIR,          
        's=f'		=>\$MIN_SCORE,   
);

#temp
$INPUT_RIS_DIR = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/VDR-BV";

my $USAGE = "\nUSAGE: $0 -m=<VDRBV_RIS_DIR> -t=<MIN_SCORE>\n" .
			"<VDRBV_RIS_DIR> ris file from PscanChip\n" .
			"(opt)<MIN_SCORE> min PscanChip score to consider (default=undef)\n";
			
unless($INPUT_RIS_DIR){
	print $USAGE;
	exit -1;
}
print "Minimum PScanChIP score set to $MIN_SCORE\n";

my $temp_ris_bed           = $INPUT_RIS_DIR . '/' . "tmp_intervals.bed";
my $temp_ris_bed_score_srt = $INPUT_RIS_DIR . '/' . "tmp_intervals_s.bed";

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



####
#2. Get all VDR-BVs from file
####
my %VDRBV_coords;
open ($instream,  q{<}, $IN_VDRBV) or die("Unable to open $IN_VDRBV : $!");
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
	$VDRBV_coords{$coord} = 1;
}
close $instream;


#delete $HASH{$key};
#############
#3. Check the vdr-bvs against the RXR:VDR ris:
#-get a bed from each risc, ordered by decreasing score
#-for each line in bed, see if vdrbv is in interval. If so, remove it from main hash and place in result hash
#############
if($RXR_VDR_RIS =~ /BVs_(.*)_(.*)_sites/){
	$motif_name = $1;
	$identifier = $2;
	$motif_string = $motif_name . '_' . $identifier;
}else{
	print STDERR "ERROR: Unable to recognise the input motif file name: $RXR_VDR_RIS. Aborting.\n";
	exit -1;
}
$motif_string = $motif_name . '_' . $identifier;
$motif_name =~ s/\-/\:\:/;
#get the length of the motif analyzed in this iteration
my $full_motif_id = $motif_name . '_' . $identifier;
my $ref = $motif{$full_motif_id};
my $MOTIF_LENGTH = scalar(@$ref); 
print STDERR "The length of the motif: $full_motif_id according to the JASPAR Pwm is $MOTIF_LENGTH\n";


open ($instream,      q{<}, $RXR_VDR_RIS) or die("Unable to open $RXR_VDR_RIS : $!");
open (my $outstream,  q{>}, $temp_ris_bed) or die("Unable to open $temp_ris_bed : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^CHR/);
	
	my ($chr,$motif_start,$motif_end,$motif_strand,$score,$site) = (split /\t/)[0,4,5,8,9,10];
	#next if (!$chr);
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
	my $bed_line = $chr . "\t" . $motif_start . "\t" . $motif_end . "\t" . $score;
	print $outstream $bed_line, "\n";
}
close $instream;
close $outstream;

#sort bed file by decreasing score
system "cat $temp_ris_bed | sort -k4,4V > $temp_ris_bed_score_srt";
#needs to be checked here

my %RESULTS;
foreach my $item (keys %VDRBV_coords){
	my $best_scoring_intersecting_pwm = check_vdrbv_intersects_pwm_interval($item, $temp_ris_bed_score_srt);
	if ($best_scoring_intersecting_pwm){
		my ($chr, $start, $stop, $score) = split("\t", $best_scoring_intersecting_pwm);
		#TODO you might want to save more than the score.
		$RESULTS{$item}{$motif_string} = $score;
		delete($VDRBV_coords{$item});
		next;
	}else{
		next;
	}
}

exit;
#do the above for all the other ris: for each ris get the best scoring intersecting pwm, if any, and then outside the ris file loop choose
#the global one with the highest score and fill results.



#reinitialise variables inside file loop
undef $identifier;
undef $motif_name;
unlink $temp_ris_bed;
unlink $temp_ris_bed_score_srt;



#subs--------------------------------------------------------------------------------------------

#############
# should return TRUE if the ris interval contains at least a VDR-BV; under otherwise
#############
sub check_vdrbv_intersects_pwm_interval{
	my ($vdrbv_coords, $filepath) = @_;
	
	my ($chr,$pos) = split("-", $vdrbv_coords);
	open (my $instream, q{<}, $filepath) or die("Unable to open $filepath : $!");
	while(<$instream>){
		chomp;
		my ($pwm_chr, $pwm_start, $pwm_end) = split("\t", $_);
		if( ($pwm_chr eq $chr) && ($pos >= $pwm_start) && ($pos <= $pwm_end) ) {
			return $_;
		}
	}
	close $instream;
	return undef;
}



#############
#4. Then check the remaining ones against all other pwms
#############
#
#chdir $INPUT_RIS_DIR;
#my @files = <Pscanchip_hg19*.ris>;
#foreach my $file (@files){
#	next if ($file eq $RXR_VDR_RIS);
#	#$file =~ s/(.*)\..*/$1/;
#	#get motif name and id from the filename
#	#eg Pscanchip_hg19_bkgGM12865_Jaspar_VDRBVs_ELK4_MA0076.1_sites.ris
#	if($file =~ /BVs_(.*)_(.*)_sites/){
#		$motif_name = $1;
#		$identifier = $2;
#	}else{
#		print STDERR "ERROR: Unable to recognise the input motif file name: $file. Aborting.\n";
#		exit -1;
#	}
#	
#	
#	
#	
#}


##heterodimers are saved by Jaspar as monomer::monomer
##I replaced the :: with a '-' in the input file name because it's not recognised by the SGE submission
##change again here:
#$motif_name =~ s/\-/\:\:/;
##get the length of the motif analyzed in this iteration
#my $full_motif_id = $motif_name . '_' . $identifier;
#my $ref = $motif{$full_motif_id};
#my $MOTIF_LENGTH = scalar(@$ref); 
#print STDERR "The length of the motif: $full_motif_id according to the JASPAR Pwm is $MOTIF_LENGTH\n";
#
#my %pwm_interval_bed;
#open ($instream,  q{<}, $input_pscanchip_ris) or die("Unable to open $input_pscanchip_ris : $!");
#while(<$instream>){
#	chomp;
#	next if($_ eq '');
#	next if($_ =~ /^CHR/);
#	
#	my ($chr,$motif_start,$motif_end,$motif_strand,$score,$site) = (split /\t/)[0,4,5,8,9,10];
#	#next if (!$chr);
#	if(  $MIN_SCORE && ($score < $MIN_SCORE) ){
#		next;
#	}
#	my $pscanchip_interval_length = ($motif_end - $motif_start);
#
#	if($pscanchip_interval_length <  $MOTIF_LENGTH){
#		$motif_end += 1;
#	}elsif($pscanchip_interval_length == $MOTIF_LENGTH){	
#		;
#	}else{
#		print STDERR "ERROR: the pscanchip motif length is LARGER than the Jaspar length. Verify.\n";
#		exit -1;	
#	}
#	my $bed_line = $chr . "\t" . $motif_start . "\t" . $motif_end;
#	$pwm_interval_bed{$bed_line} = 1;
#}
#close $instream;
#
#open (my $outstream,  q{>}, $temp_pwm_bed) or die("Unable to open $temp_pwm_bed : $!");
#foreach my $item (keys %pwm_interval_bed){ print $outstream $item, "\n"; }
#close $outstream;
#
##sort bed
#system "cat $temp_pwm_bed | sort -k1,1V -k2,2n | uniq > $temp_pwm_bed_sorted";
##get closeness data
#system "$BEDTOOLS closest -D \"ref\" -a $input_variants -b $temp_pwm_bed_sorted -g $CHROMSIZES > $data_closest";
#unlink $input_variants if($INTERVALS);
##from here I want two files
##one, pure counts, to build histograms in R
##two, x,y pairs, to build point plots in R
#
##get pure counts
#system "cat $data_closest  | cut -f 7 > $data_histogram";
##get x,y pairs, where x=distance and y=count at that distance
#system "cat $data_closest  | cut -f 7 | sort -k1,1n | uniq -c | sed \'s/^ *//g;s/ /\t/\'  > $data_counts";
#
#unlink $temp_pwm_bed;
#unlink $temp_pwm_bed_sorted;
#unlink $data_closest;





