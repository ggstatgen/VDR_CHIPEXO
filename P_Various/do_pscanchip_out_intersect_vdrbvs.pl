#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#/11/12/2015
#I FOUND OUT THAT THE INTERVAL IN THE PSCANCHIP ris FILE CAN BE WRONG, DEPENDING ON THE TF considered!!!
#the results of this file are probably redundant.

#Actually, this is probably redundant. Funseq2 will already intersect the .ris intervals with the VDR-Bv, so you don't need this PRIOR to funseq2
#This is useful to have an initial approximation of how many interesections you get AT DIFFERENT SCORES

#30/11/2015
#REVIEW point: RXR:VDR motifs under VDR-BVs

#Given the results from PScanChIP starting from a bed of VDR-BVs, I need to know exactly which returned intervals overlap a RXR:VDR motif.
#I cannot say this from the PscanChIP output, because the program extends ever 1nt VDR-BV left and right by 150nt.

#This program will return a subset of the initial file having only the rows where the VDR-BV is INSIDE the motif found.

#NOTE 




#format is as follows
#CHR     REG_START       REG_END REG_STRAND      ABS_SITE_START  ABS_SITE_END    REL_SITE_START  REL_SITE_END    SITE_STRAND     SCORE   SITE
#chr9    131645230       131645379       +       131645239       131645253       -66     -53     +       0.989252        GGGTCATGGAGTTCA
#chr6    34302932        34303081        +       34302996        34303010        -11     2       -       0.989252        TGAACTCCATGACCC

#ris coordinates are -1 of the ones needed by UCSC. 
#If you want the UCSC BROWSER ones, add 1 to both
#If you want a standard BED, add 1 to end
#If you want the BIGWIGTOWIG, add 1 to end

my $infile_vdrbvs;
my $infile_pscanchip;
my $MIN_SCORE;
GetOptions(
        'v=s'        =>\$infile_vdrbvs,
        'p=s'        =>\$infile_pscanchip, 
        'm=f'        =>\$MIN_SCORE
);

#$infile_vdrbvs = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/SUPPL_DATA_Output_noDBRECUR_REP.vcf";
#$infile_pscanchip = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/PSCANCHIP_motifs/VDR-rBV/Pscanchip_hg19_bkgGM12865_Jaspar_VDRrBVs_NR1h3RXRA_sites.ris";

if(!$infile_vdrbvs){
	print "USAGE: do_pscanchip_out_intersect_vdrbvs.pl -v=<INFILE_VDRBVs> -p=<INFILE_PSCANCHIP> (opt)-m=<MINSCORE>\n";
    print "<INFILE_VDRBVs> .vcf file out of funseq, with the VDR-BVs\n";
    print "<INFILE_PSCANCHIP> .res/ris file with VDR:RXR intervals out of Pscanchip [hg19]";
    print "OPTIONAL - <MINSCORE> lower threshold on score (eg 0.8) (default: none)\n";
    exit 1;
}
if(!$infile_pscanchip){
    print "USAGE: do_pscanchip_out_intersect_vdrbvs.pl -v=<INFILE_VDRBVs> -p=<INFILE_PSCANCHIP> (opt)-m=<MINSCORE>\n";
    print "<INFILE_VDRBVs> .vcf file out of funseq, with the VDR-BVs\n";
    print "<INFILE_PSCANCHIP> .res/ris file with VDR:RXR intervals out of Pscanchip [hg19]";
    print "OPTIONAL - <MINSCORE> lower threshold on score (eg 0.8) (default: none)\n";
    exit 1;
}
print STDERR "THRESHOLDING ON SCORE: $MIN_SCORE\n" if($MIN_SCORE);

#save the vdr chr-pos in a data structure
my %variant_coords;
open (my $instream,  q{<}, $infile_vdrbvs) or die("Unable to open $infile_vdrbvs : $!");
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


#analyze the intervals
#each VDR-BV MUST appear in at least one [REG_START, REG_END] interval
#so check this interval against the VDRs above
#1) map vdr-bv to VDR-BV-interval
#ONCE you found the VDR-BV-interval, check if the VDR-BV is in the [ABS_SITE_START,ABS_SITE_END]

my %dataline; #save line in here, print later
open ($instream,  q{<}, $infile_pscanchip) or die("Unable to open $infile_pscanchip : $!");
#print track
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^CHR/);
	
	$dataline{$_} = 1;
}
close $instream;

#do the comparisons
my %vdrbv_seen_in_the_intervals;
my %final_data;

foreach my $item (keys %variant_coords){
	my ($vdrbv_chr, $vdrbv_pos) = split('-', $item);

	foreach my $vdrbv_interval (keys %dataline){
		my ($pscanchip_chr, 
			$reg_start, $reg_end, 
			$abs_site_start, $abs_site_end, 
			$rel_site_start, $rel_site_end, $score) = (split /\t/, $vdrbv_interval)[0,1,2,4,5,6,7,9];	
		
		#turn into bed format
		$reg_end += 1;
		$abs_site_end += 1;
		#what to do with the relative position? leave them as they are for now

		next if (!$pscanchip_chr);
		next if(  $MIN_SCORE && ($score < $MIN_SCORE) );
		next unless ($pscanchip_chr eq $vdrbv_chr);
		next unless ( ($vdrbv_pos >= $reg_start) && ($vdrbv_pos <= $reg_end )  );
		next unless ( ($vdrbv_pos >= $abs_site_start) && ($vdrbv_pos <= $abs_site_end )  );
		
		$vdrbv_seen_in_the_intervals{$item} = 1;
		$final_data{$vdrbv_interval} = 1;
	}
}

#print STDOUT "#CHR\tREG_START\tREG_END\tREG_STRAND\tABS_SITE_START\tABS_SITE_END\tREL_SITE_START\tREL_SITE_END\tSITE_STRAND\tSCORE\tSITE\n";
foreach my $line (keys %final_data){
	print STDOUT $line, "\n";
}




