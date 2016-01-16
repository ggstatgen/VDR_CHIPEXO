#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#(old)
###FIGURE 5D
#Script to map VDR-BVs and BACKGROUND SNPs to D' LD blocks

#derived from do_funseq_GAT_LDclump_FG_and_BG.pl
#this is a conceptually much simpler version replacing all variants with the intervals they're in

#for example:
#    1                2  3         4
#----.----------------.--.---------.--------
#---_____-----------________-----------______
#result (output bed interval):
#---_____-----------________-------.--------

#if a variant is not in an LD block, report the variant's interval
#else, report the LD block

my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";
my $INPUT_LD = "/net/isi-mirror/1000_Genomes/POPULATION_PANELS/PLINK_HAPLOTYPE_BLOCKS/plink_CEU_LD_b37_noindels_MAF0.01.bed.gz"; #general LD blocks #b37
my $FOREGROUND_SET;
my $in_foreground; my $in_background;
GetOptions(
        'fg=s'      =>\$in_foreground,
        'bg=s'      =>\$in_background,
        'fgs=s'     =>\$FOREGROUND_SET
);
#CONVERT ALL TO HG19
if(!$in_foreground){
	print "USAGE: do_funseq_GAT_LDclump_FG_and_BG_2BLOCKS.pl -fg=<FG> -bg=<BG> -fgs=<ALL|HCBV|REP|ENH>\n";
    print "<FG> vcf output of funseq2 with noRECUR\n";
    print "<BG> vcf.gz file containing your choice of variants for GAT background\n";
    print "fgs: foreground set - one of ALL (all VDR-BVs), HCBV (RECUR), REP (only replicating), ENH (in enhancer)\n";
    exit 1;
}
if(!$in_background){
	print "USAGE: do_funseq_GAT_LDclump_FG_and_BG_2BLOCKS.pl -fg=<FG> -bg=<BG> -fgs=<ALL|HCBV|REP|ENH>\n";
    print "<FG> vcf output of funseq2 with noRECUR\n";
    print "<BG> vcf.gz file containing your choice of variants for GAT background\n";
    print "fgs: foreground set - one of ALL (all VDR-BVs), HCBV (RECUR), REP (only replicating), ENH (in enhancer)\n";
    exit 1;
}
if(!$FOREGROUND_SET){
	print "USAGE: do_funseq_GAT_LDclump_FG_and_BG_2BLOCKS.pl -fg=<FG> -bg=<BG> -fgs=<ALL|HCBV|REP|ENH>\n";
    print "<FG> vcf output of funseq2 with noRECUR\n";
    print "<BG> vcf.gz file containing your choice of variants for GAT background\n";
    print "fgs: foreground set - one of ALL (all VDR-BVs), HCBV (RECUR), REP (only replicating), ENH (in enhancer)\n";
    exit 1;
}
unless($FOREGROUND_SET eq 'ALL' || $FOREGROUND_SET eq 'HCBV' || $FOREGROUND_SET eq 'REP' || $FOREGROUND_SET eq 'ENH'){
	print "USAGE: do_funseq_GAT_LDclump_FG_and_BG_2BLOCKS.pl -fg=<FG> -bg=<BG> -fgs=<ALL|HCBV|REP|ENH>\n";
    print "<FG> vcf output of funseq2 with noRECUR\n";
    print "<BG> vcf.gz file containing your choice of variants for GAT background\n";
    print "fgs: foreground set - one of ALL (all VDR-BVs), HCBV (RECUR), REP (only replicating), ENH (in enhancer)\n";
    exit 1;
}
my($basename, $directory) = fileparse($in_foreground);
$basename =~ s/(.*)\..*/$1/;

########
#outputs
########
my $out_temp_fg_bed    = $directory . $basename . '_temp_FG.bed'; 
my $out_temp_bg_bed    = $directory . $basename . '_temp_BG.bed';

my $out_temp_fg_ol     = $directory . $basename . '_temp_FG_ol.bed';
my $out_temp_fg_nool   = $directory . $basename . '_temp_FG_nool.bed';
my $out_temp_fg_nool_pr = $directory . $basename . '_temp_FG_nool_proc.bed';
     
my $out_temp_bg_ol     = $directory . $basename . '_temp_BG_ol.bed';
my $out_temp_bg_nool   = $directory . $basename . '_temp_BG_nool.bed';
my $out_temp_bg_nool_pr = $directory . $basename . '_temp_BG_nool_proc.bed';

my $output_fg                =  $directory . 'LDCLUMP_BLOCK_FOREGROUND_' . $FOREGROUND_SET . '_' . $basename . '.bed.gz';
my $output_bg                =  $directory . 'LDCLUMP_BLOCK_BACKGROUND_'                         . $basename . '.bed.gz';

my %data_fg;
my %data_bg;
my %no_overlap_fg;
my %no_overlap_bg;
################
#1.1 get a bed of the variants in the foreground set from the output of funseq
################
open (my $instream,   q{<}, $in_foreground) or die("Unable to open $in_foreground : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');	
	next if($_ =~ /^\#/);
	my $info_ncenc; my $info_motifbr; my $info_recur;
	
	my ($chr, $pos, $info) = (split /\t/)[0,1,7]; #chr is hg19	

	my @info = split(";", $info);
	foreach my $item (@info){
		$info_ncenc = $item if($item =~ /^NCENC/); 
		$info_motifbr = $item if($item =~ /^MOTIFBR/); #any motif
		$info_recur = $item  if($item =~ /^RECUR/);			
	}
	$chr =~ s/chr(.*)/$1/; #this needs to be b37
	my $vdrhcbv = $chr . '-' . $pos;
	
	#FLAGS
	if($FOREGROUND_SET eq 'ALL'){
		$data_fg{$vdrhcbv} = 1;
		next;
	}elsif($FOREGROUND_SET eq 'REP'){
		if(!$data_fg{$vdrhcbv}){
			$data_fg{$vdrhcbv} = 1;
		}else{
			$data_fg{$vdrhcbv} += 1;
		}
		next;	
	}elsif($FOREGROUND_SET eq 'HCBV'){
		if($info_recur){
			$data_fg{$vdrhcbv} = 1;			
		}else{
			next;
		}
	}elsif($FOREGROUND_SET eq 'ENH'){
		next unless($info_ncenc);
		if($info_ncenc =~ /Enhancer/){
			$data_fg{$vdrhcbv} = 1;				
		}else{
			next;				
		}
	}else{
		print STDERR "Error: unrecognised: $FOREGROUND_SET. Aborting..\n";
		exit -1;
	}
}
close $instream;
	
open (my $outstream,  q{>}, $out_temp_fg_bed) or die("Unable to open $out_temp_fg_bed: $!");		
foreach my $item (sort keys %data_fg){
	
	if($FOREGROUND_SET eq 'REP'){
		next unless($data_fg{$item} > 1);
	}
	my ($chr, $pos) = split("-", $item);
	print $outstream $chr . "\t" . ($pos - 1) . "\t" . $pos . "\n";
}
close $outstream;

################
#1.2 get a bed with the background variants. Remove indels
################
tie *FILE,   'IO::Zlib', $in_background, "rb";
while (<FILE>)	{ 
	next if($_ eq '');	
	next if($_ =~ /^\#/);
	chomp;
	my ($chr, $pos, $id, $info) = (split /\t/)[0,1,2,7];	
	
	#remove cases which have "esv" instead of rs id
	next if($id =~ /^esv/);
	#remove INDELS
	next if($info =~ /INDEL/);
	my $variant_bkg = $chr . '-' . $pos; 
	$data_bg{$variant_bkg} = 1; 
}
close FILE;

open ($outstream,  q{>}, $out_temp_bg_bed) or die("Unable to open $out_temp_bg_bed: $!");		
foreach my $item (sort keys %data_bg){
	my ($chr, $pos) = split("-", $item);
	print $outstream $chr . "\t" . ($pos - 1) . "\t" . $pos . "\n";
}
close $outstream;

#variants which are in LD blocks
system "$BEDTOOLS intersect -wo -a $INPUT_LD -b $out_temp_fg_bed > $out_temp_fg_ol"; #fg
system "$BEDTOOLS intersect -wo -a $INPUT_LD -b $out_temp_bg_bed > $out_temp_bg_ol"; #bg
#variants which are NOT in LD blocks
system "$BEDTOOLS intersect -v -a $out_temp_fg_bed -b $INPUT_LD  > $out_temp_fg_nool"; #fg
system "$BEDTOOLS intersect -v -a $out_temp_bg_bed -b $INPUT_LD  > $out_temp_bg_nool"; #bg
#convert the latter in the format of the first and cat, sort, zip

#FG--------------
convert_nool_2_ol($out_temp_fg_nool, \%no_overlap_fg);
open ($outstream,  q{>}, $out_temp_fg_nool_pr) or die("Unable to open $out_temp_fg_nool_pr : $!");
foreach my $item (keys %no_overlap_fg){ print $outstream $item, "\n"; }
close $outstream;
system "cat $out_temp_fg_ol $out_temp_fg_nool_pr | sort -k1,1V -k2,2g | gzip -c > $output_fg";

#BG--------------
convert_nool_2_ol($out_temp_bg_nool, \%no_overlap_bg);
open ($outstream,  q{>}, $out_temp_bg_nool_pr) or die("Unable to open $out_temp_bg_nool_pr : $!");
foreach my $item (keys %no_overlap_bg){ print $outstream $item, "\n"; }
close $outstream;
system "cat $out_temp_bg_ol $out_temp_bg_nool_pr | sort -k1,1V -k2,2g | gzip -c > $output_bg";

unlink $out_temp_fg_bed; 
unlink $out_temp_bg_bed;
unlink $out_temp_fg_ol;
unlink $out_temp_fg_nool;
unlink $out_temp_fg_nool_pr;
unlink $out_temp_bg_ol;
unlink $out_temp_bg_nool;
unlink $out_temp_bg_nool_pr;

#---------------------subs--------------------------
sub convert_nool_2_ol{
	my ($file, $hash) = @_;
	
	open (my $instream,  q{<}, $file) or die("Unable to open $file: $!");
	while(<$instream>){
		chomp;
		my @fields = split("\t", $_);
		my $new_bed_line = $fields[0] . "\t" . ( $fields[1]- 1 ) . "\t" . $fields[1] . "\t(1)" . $fields[2] . "\t" .  '-';
		my $new_line = $new_bed_line . "\t" . $_ . "\t" . '0';
		$$hash{$new_line} = 1;
	}
	close $instream;
	return;
}