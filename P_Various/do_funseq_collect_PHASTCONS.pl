#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#20/1/2015
#Figure 3 show low numbers of variants on two highly conserved nucleotides.
#Chris asked to see the phastcons scores for the motif sequences.
#This should get phastcons scores and return R data so that I can plot some box plots and hopefully show that these nucleotides are strongly conserved

#This phastcons data is in bigwig format, convert as follows
#bigWigToWig - Convert bigWig to wig.  This will keep more of the same structure of the
#original wig than bigWigToBedGraph does, but still will break up large stepped sections
#into smaller ones.
#usage:
#   bigWigToWig in.bigWig out.wig
#options:
#   -chrom=chr1 - if set restrict output to given chromosome
#   -start=N - if set, restrict output to only that over start
#   -end=N - if set, restict output to only that under end
#   -udcDir=/dir/to/cache - place to put cache for remote bigBed/bigWigs

my $INFILE_RIS; 
my $regions_bed;
#my $PHASTCONS_BW = "/net/isi-mirror/ucsc/hg19/primatesPhast.bw";
my $PHASTCONS_BW = "/net/isi-mirror/ucsc/hg19/vertebratesPhast.bw";
my $UCSC_bigwig2wig = "/net/isi-scratch/giuseppe/tools/UCSC_tools/bigWigToWig";
#This contains strand information, which I need to obtain the start/stop positions to query bigWigToWig with
GetOptions(
        'i=s'       =>\$INFILE_RIS,
        'regions=s' =>\$regions_bed
); 
#$INFILE_RIS = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_motif_scan_pscanchip/Pscanchip_occurrences_VDRRXR.all.ris";
#$INFILE_RIS = "/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_motif_scan_pscanchip/temp.ris";
if(!$INFILE_RIS){
	print "USAGE: do_funseq_collect_PHASTCONS.pl -i=<INFILE> -regions=<REGIONSET>\n";
	print "<INFILE> pscanchip .ris interval file with MOTIF INTERVALS and strand info\n";
    print "<REGIONSET> bed file with VDR:RXR PEAK intervals to restrict the search to\n";
    print "NOTE: EDIT script to change phastcons data [vertebrate|primate]\n";
    exit 1;
}
if(!$regions_bed){
	print "USAGE: do_funseq_collect_PHASTCONS.pl -i=<INFILE> -regions=<REGIONSET>\n";
	print "<INFILE> pscanchip .ris interval file with MOTIF INTERVALS and strand info\n";
    print "<REGIONSET> bed file with VDR:RXR PEAK intervals to restrict the search to\n";
    print "NOTE: EDIT script to change phastcons data [vertebrate|primate]\n";
    exit 1;
}
my($basename, $directory) = fileparse($regions_bed);
$basename =~ s/(.*)\..*/$1/;
my $temp_wig = $directory . $basename . '_temp.wig';
my $output  = $directory  . $basename . '_phastcons.Rdata';

################
#0 bed file with peaks (eg class I peaks, enhancer peaks) save their coordinates in hash
################
my %regions_bed;
open (my $instream,  q{<}, $regions_bed) or die("Unable to open $regions_bed : $!");
while(<$instream>){
	chomp;
	my ($chr, $start, $stop) = (split /\t/)[0,1,2];
	my $interval = $start . '-' . $stop;
	$regions_bed{$chr}{$interval} = 1;
}
close $instream;

#####################
#collect unique motif intervals; reformat intervals to agree with bwTOwig
#####################
my %motif_to_strand;
open ($instream,  q{<}, $INFILE_RIS) or die("Unable to open $INFILE_RIS : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^CHR/);
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
	my $bw2w_start = $fields[4]; #bw2w will count from this + 1
	my $bw2w_stop = $fields[5] + 1;
	
	my $identifier = $fields[0] . '-' . $bw2w_start . '-' . $bw2w_stop;     #chr7-5013558-5013573
	$motif_to_strand{$identifier} = $fields[8];
}
close $instream;

####
#go through the collected motif intervals and get phastcons scores
####
my %score_matrix; #each key is a motif_id (chr-start-end) and each value is a 15-ple of phastcons scores, from smaller to bigger coordinate, split by tab
foreach my $item (sort keys %motif_to_strand){
	my @phastcons_scores;
	my ($chr, $start, $stop) = split("-", $item);
	
	#this motif's coordinates have to fit into a corresponding peak coordinate - othewise, skip motif
	next unless(defined check_coords_in_bed($chr, $start, $stop));
	
	system "$UCSC_bigwig2wig $PHASTCONS_BW $temp_wig -chrom=$chr -start=$start -end=$stop";
	#collect data from wig file, save in hash/array, unlink wig file, next
	#position will be inverted if strand = -
	if ( -z $temp_wig) {
		print STDERR "Warning: no Phastcons data for interval $item. Skipping..\n";
		next;
	}
	#wig output MUST have 16 columns (1 header + 15 values)
	#in some outputs there is more than one header
	my $lines = `cat $temp_wig | grep -v "fixedStep" | wc -l`;
	unless($lines =~ /^15/){
		print STDERR "ERROR: motif interval size: $lines is not 15. Chromosome: $chr; Start: $start; end: $stop. Skipping..\n";
		system "cat $temp_wig";
		#exit -1;
		next
	}
	#wig output MUST look like this:
#	fixedStep chrom=chr9 start=131645239 step=1 span=1
#	0.047
#	..x15
	#TODO you need to deal with gaps!!
	open (my $tempstream,  q{<}, $temp_wig) or die("Unable to open $temp_wig : $!");
	while(<$tempstream>){
		chomp;
		next if($_ =~ /^fixedStep/);
		if($_ eq ''){
			push(@phastcons_scores,'NA');
		}else{
			push(@phastcons_scores,$_);
		}
	}
	@phastcons_scores = reverse(@phastcons_scores) if($motif_to_strand{$item} eq '-');
	$score_matrix{$item} = join("\t",@phastcons_scores);
	close $tempstream;
	unlink $temp_wig;
}
open (my $outstream,  q{>}, $output) or die("Unable to open $output : $!");
#print $outstream "CHR\tSTART\tSTOP\tp1\tp2\tp3\tp4\tp5\tp6\tp7\tp8\tp9\tp10\tp11\tp12\tp13\tp14\tp15\n";
print $outstream "ID\tp1\tp2\tp3\tp4\tp5\tp6\tp7\tp8\tp9\tp10\tp11\tp12\tp13\tp14\tp15\n";
foreach my $item (sort keys %score_matrix){
	#print $outstream $chr . "\t" . $start . "\t" . $stop . "\t" . $score_matrix{$item} . "\n";
	print $outstream  $item . "\t" .  $score_matrix{$item} . "\n";
}
close $outstream;


###########
#subs
###########
sub check_coords_in_bed{
	my ($mtf_chr, $mtf_start, $mtf_stop) = @_;
	
	foreach my $peak_interval (keys %{ $regions_bed{$mtf_chr} }){
		my ($peak_start, $peak_end) = split('-', $peak_interval);
		if( ($mtf_stop <= $peak_end)  &&  ($mtf_start >= $peak_start) ){
			return 1;
		}
	}
	return undef;
}