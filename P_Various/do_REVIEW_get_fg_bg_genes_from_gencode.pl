#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#24/1/2016
#I want to see if the VDR-BVs tend to be in the vicinity of genes known to be VITD differentially regulated, taken from Ramagopalan 2010:
#http://genome.cshlp.org/content/suppl/2010/08/24/gr.107920.110.DC1/Supplementary_Table_2.xls
#These are b36, my VDR-BVs are hg19
#I tried to lift over, I could not find a perfect correspondence in the coordinates by looking at Ensembl 
#http://grch37.ensembl.org/index.html


#So went to gencode and downloaded the last b37 version, gencode 19.
#This scripts will extract:
#1) a background of all protein coding genes (GAT background)
#2) the foreground of genes from the table above

#INPUTS
#1 list of Ensembl ID from the table in the vit d paper
#protein coding genes from gencode

#OUTPUTS
#2 beds in hg19, one for the foreground, the other for the background

#Then you will extend your VDR-BVs by 1MB and do a GAT

my $input_gencode;
my $input_list;
my $lifted_over_in;

my $VCFSORT = `which vcf-sort`; chomp $VCFSORT;
my $GZIP = `which gzip`; chomp $GZIP;

GetOptions(
        'i=s'      =>\$input_list,
		'g=s'      =>\$input_gencode,
		'l=s'      =>\$lifted_over_in
); 
my $USAGE = "\nUSAGE: do_REVIEW_get_fg_bg_genes_from_gencode.pl -i=<INFILE_LIST> -g=<GENCODE_PC_GENES> (opt)-l=<LIFTED_OVER_BED>\n" .
			"<INFILE_LIST> list of foreground ENSID to convert extracted from VITD Chipseq supplementary\n" .
			"<GENCODE_PC_GENES> gtf gz file with GENCODE v19 entries of the type 'gene' and 'protein_coding'\n" . 
			"(optional)<LIFTED_OVER_BED> hg19 bed file lifted over from the hg18 bed - used to get coords if ENS ID not found in gencode.\n";
			
unless($input_list && $input_gencode){
	print STDERR $USAGE, "\n";
	exit 1;
}
my($basename, $directory) = fileparse($input_list);
my($gen_basename, $gen_directory) = fileparse($input_gencode);
$basename =~ s/(.*)\..*/$1/;
$gen_basename =~ s/(.*)\..*/$1/;
my $output_fg = $directory . $basename . '_fg_hg19.bed';
my $output_bg = $directory . $basename . '_bg_hg19.bed';
my $output_not_found = $directory. $basename . '_unmapped_to_' . $gen_basename . '.txt';


#create background file
#input:

#chr1
#HAVANA
#gene
#69091
#70008
#.
#+
#.
#gene_id "ENSG00000186092.4"; transcript_id "ENSG00000186092.4"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F5"; level 2; havana_gene "OTTHUMG00000001094.1";

my %lifted_over_map;
if($lifted_over_in){
	open (my $instream,  q{<}, $lifted_over_in) or die("Unable to open $lifted_over_in : $!");
	while(<$instream>){
		chomp;
		my @fields = split("\t", $_);
		print $fields[3], "\n";
		my ($gene_name, $ens_id) = split("|", $fields[3]);
		$lifted_over_map{$ens_id}{CHR} = $fields[0];
		$lifted_over_map{$ens_id}{START}= $fields[1];
		$lifted_over_map{$ens_id}{STOP}= $fields[2];
		$lifted_over_map{$ens_id}{GENE}= $gene_name;	
	}
	close $instream;
}



my %fgdata;
open (my $instream,  q{<}, $input_list) or die("Unable to open $input_list : $!");
open (my $outstream,  q{>}, $output_not_found) or die("Unable to open $output_not_found : $!");
while(<$instream>){
	chomp $_;
	my $ens_id = $_;
	my $grep_search = "zcat $input_gencode | grep \"$ens_id\"";
	my $line = `$grep_search`;
	
    if(!$line || $line eq '' ){
    	print STDERR "Warning: grepping $ens_id against the GENCODE data found nothing:";
    	if($lifted_over_in){
    		if($lifted_over_map{$ens_id}){
    			print STDERR "Found $ens_id in lifted over file..\n";
    			my $bed_line = $lifted_over_map{$ens_id}{CHR}    . "\t" .
    							$lifted_over_map{$ens_id}{START} . "\t" .
    							$lifted_over_map{$ens_id}{STOP}  . "\t" .
    							'orig_coords_' . $ens_id . '-' . $lifted_over_map{$ens_id}{GENE};
    			$fgdata{$bed_line} = 1;
    			next;
    		}else{
    			print STDERR "Not found $ens_id in lifted over file..\n";
    			print $outstream $ens_id, "\n";
    			next;
    		}
    	}else{
    		print STDERR "saving this in dump file\n";
    		print $outstream $ens_id, "\n";
    		next;
    	}
    }
    my @lines = split("\n", $line);
    if(scalar @lines > 1){
    	print STDERR "Warning: more than one line retrieved for $ens_id:\n";
    	foreach my $item (@lines){ print $item, "\n"; }
    	exit -1;
    }
    my @fields = split("\t", $lines[0]);
    next unless($fields[3]);
	next unless($fields[4]);
 	if($fields[2] ne 'gene'){
		print STDERR "Warning: this line does not contain a gene: $lines[0]. Skipping..\n";
		next;
	}

	my @inner_fields = split(';', $fields[8]);
	my $bed_name = $ens_id . ';' . $inner_fields[4];
	#ADJUST COORDS 
	my $bed_line = $fields[0] . "\t" . ($fields[3] - 1) . "\t" . $fields[4] . "\t" . $bed_name; 
  	$fgdata{$bed_line} = 1;
}
close $instream;
close $outstream;

#print foreground
open ($outstream,  q{>}, $output_fg) or die("Unable to open $output_fg : $!");
foreach my $item (keys %fgdata){ print $outstream $item, "\n"; }
close $outstream;

my %bgdata;
tie *FILE,   'IO::Zlib', $input_gencode, "rb";
while (<FILE>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^\#/);
	next if($_ =~ /^\#/);
	
	my @fields = split("\t", $_);
	#general check
	if($fields[2] ne 'gene'){
		print STDERR "Warning: this line does not contain a gene: $_. Skipping..\n";
		next;
	}
	next unless($fields[3]);
	next unless($fields[4]);

	my @inner_fields = split(';', $fields[8]);
	my $bed_name = $inner_fields[0] . '; ' . $inner_fields[4];
	my $bed_line = $fields[0] . "\t" . ($fields[3] - 1) . "\t" . $fields[4] . "\t" . $bed_name; 
  	$bgdata{$bed_line} = 1;
}
close FILE;

#print background
open ($outstream,  q{>}, $output_bg) or die("Unable to open $output_bg : $!");
foreach my $item (keys %bgdata){ print $outstream $item, "\n"; }
close $outstream;


