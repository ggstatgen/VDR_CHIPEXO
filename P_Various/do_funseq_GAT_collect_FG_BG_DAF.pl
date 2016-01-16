#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Zlib;
use List::MoreUtils qw(:all);

#derived from do_funseq_RAM_collect_DAF_GAT_background.pl

#given a bed of foreground vars and a bed of background vars, maybe from do_funseq_GAT_LDclump_FG_andBG.pl, create some input for R to determine if they're frequency matched

my $in_foreground; 
my $in_background;
my $INPUT_1KG = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/d_GAT_BACKGROUNDS/GAT_BACKGROUND_ALL.vcf.gz";
my $FILE_TYPE;
GetOptions(
        'fg=s'      =>\$in_foreground,
        'bg=s'      =>\$in_background,
        't=s'       =>\$FILE_TYPE
);
#CONVERT ALL TO HG19
if(!$in_foreground){
	print "USAGE: do_funseq_GAT_collect_FG_BG_DAF.pl -fg=<FG> -bg=<BG> -t=<bed|vcf>\n";
    print "<FG> b37 gzipped file with foreground coordinates\n";
    print "<BG> b37 gzipped file with background coordinates\n";
    print "<t> whether FG and BG are bed.gz or vcf.gz files\n";
    exit 1;
}
if(!$in_background){
	print "USAGE: do_funseq_GAT_collect_FG_BG_DAF.pl -fg=<FG> -bg=<BG> -t=<bed|vcf>\n";
    print "<FG> b37 gzipped file with foreground coordinates\n";
    print "<BG> b37 gzipped file with background coordinates\n";
    print "<t> whether FG and BG are bed.gz or vcf.gz files\n";
    exit 1;
}
if(!$FILE_TYPE){
	print "USAGE: do_funseq_GAT_collect_FG_BG_DAF.pl -fg=<FG> -bg=<BG> -t=<bed|vcf>\n";
    print "<FG> b37 gzipped file with foreground coordinates\n";
    print "<BG> b37 gzipped file with background coordinates\n";
    print "<t> whether FG and BG are bed.gz or vcf.gz files\n";
    exit 1;	
}
unless($FILE_TYPE eq 'bed' || $FILE_TYPE eq 'vcf'){
	print "USAGE: do_funseq_GAT_collect_FG_BG_DAF.pl -fg=<FG> -bg=<BG> -t=<bed|vcf>\n";
    print "<FG> b37 gzipped file with foreground coordinates\n";
    print "<BG> b37 gzipped file with background coordinates\n";
    print "<t> whether FG and BG are bed.gz or vcf.gz files\n";
    exit 1;		
}


my @variants_fg_eur; 
my @variants_fg_afr;
my @variants_fg_glob;
my @variants_bg_eur; 
my @variants_bg_afr;
my @variants_bg_glob;

my %variants_1kg;
#####################
#1 -build  hash to map chr-pos to DAF-CEU, DAF-YRI, DAF-GLOB
#####################
#before saving the frequencies, you need to to KNOW if the ancestral allele is the ref or the alt
#if the ancestral allele is the ref, save the frequency as is
#if the ancestral allele is the alt, the frequency you have is for the ancestral. The derived will be (1 - freq)
tie *FILE,   'IO::Zlib', $INPUT_1KG, "rb";
while (<FILE>)	{ 
	chomp;
	next if($_ eq '');
	next if($_ =~ /^#/);
	my $ALTF_AFR; my $ALTF_EUR; my $ALTF_GLOB;
	my $FLAG; # set to one if the ancestral is the alternate	

	my @fields = split("\t", $_);
	#skip indels
	next if(length($fields[3]) > 1);
	next if(length($fields[4]) > 1);
	my $key = $fields[0] . '-' . $fields[1];
	#get derived allele frequencies
	my $ref = $fields[3]; my $alt = $fields[4]; my $info = $fields[7];
	
	#get ancestral allele info=====================
	my @info = split(";", $fields[7]);
	if($info[0] =~ /^AA=(.*)/){
		my $anc = uc($1);
		if($anc eq $alt){
			$FLAG = 1;		
		}elsif($anc eq $ref){			
		}elsif( ($anc eq '') or ($anc eq 'N') or ($anc eq '.') or ($anc eq '-') ){
			next;
		}else{ 
			next;			
		}
	}else{
		next;		
	}
	#get allele frequencies========================
	foreach my $item (@info){
		if($item =~ /^AFR_AF=([0-9]+\.[0-9]+)/){
			$ALTF_AFR = $1;				
		}	
		if($item =~ /^EUR_AF=([0-9]+\.[0-9]+)/){
			$ALTF_EUR = $1;	
		}
		if($item =~ /^AF=([0-9]+\.[0-9]+)/){
			$ALTF_GLOB = $1;
		}
	}		
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
	if($ALTF_GLOB){
		if($FLAG){
			$variants_1kg{$key}{DAFGLOB} = (1 - $ALTF_GLOB);
		}else{
			$variants_1kg{$key}{DAFGLOB} = $ALTF_GLOB;
		}
	}
}
close FILE;

#you now have a hash with frequencies for CEU,YRI and global. Go through foreground and background and collect frequencies

###############
#2 get DAF for the three ethnicity groups
###############
get_DAF($in_foreground,\@variants_fg_eur,\@variants_fg_afr,\@variants_fg_glob);
get_DAF($in_background,\@variants_bg_eur,\@variants_bg_afr,\@variants_bg_glob);


print "ETHNICITY\tDAF_FG\tDAF_BG\n";
#create tsv
#excellent solution found here:
#http://www.perlmonks.org/?node_id=446821
my $result_ceu     = each_array(@variants_fg_eur,  @variants_bg_eur);
my $result_yri     = each_array(@variants_fg_afr,  @variants_bg_afr);
my $result_glob    = each_array(@variants_fg_glob, @variants_bg_glob);

while ( my ($a, $b) = $result_ceu->() ){
	$a = 'NA' unless($a);
	$b = 'NA' unless($b);
    print "CEU\t$a\t$b\n";
}
while ( my ($a, $b) = $result_yri->() ){
	$a = 'NA' unless($a);
    $b = 'NA' unless($b);
    print "YRI\t$a\t$b\n";
}
while ( my ($a, $b, $c, $d) = $result_glob->() ){
	$a = 'NA' unless($a);
	$b = 'NA' unless($b);
    print "GLOB\t$a\t$b\n";
}


#subs------------------------------------------------------------

sub get_DAF{
	my ($file, $array_eur, $array_afr, $array_glob) = @_;
	
	#map file to hash first to avoid any duplicate positions
	my %inhash;
	tie *FILE,   'IO::Zlib', $file, "rb";
	while (<FILE>)	{ 
		chomp;
		next if($_ eq '');
		
		my $chr; my $pos;
		if($FILE_TYPE eq 'vcf'){
			($chr, $pos) = (split /\t/)[0,1];			
		}elsif($FILE_TYPE eq 'bed'){
			($chr, $pos) = (split /\t/)[0,2];
		}else{
			print STDERR "ERROR: file type: $FILE_TYPE not recognised\n";
			exit -1;
		}
		my $key = $chr . '-' . $pos;
		$inhash{$key} = 1;
	}
	close FILE;
	
	foreach my $item (sort keys %inhash){
		push(@$array_eur, $variants_1kg{$item}{DAFEUR})  if($variants_1kg{$item}{DAFEUR});
		push(@$array_afr, $variants_1kg{$item}{DAFAFR})  if($variants_1kg{$item}{DAFAFR});
		push(@$array_glob,$variants_1kg{$item}{DAFGLOB}) if($variants_1kg{$item}{DAFGLOB});		
	}
	return;	
}