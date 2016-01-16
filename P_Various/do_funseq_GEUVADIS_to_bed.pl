#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#Chris wants me to cross reference VDR-BVs and eQTLs from the geuvadis data:
#I think in the abstract we can say that we predict XX genes to be differentially under the control of VDR in the human population (i.e. VDR-BVs coincident with LCL #eQTLs, despite being calcitriol-minus) and YY genes that are both under the control of VDR and contribute to disease susceptibility (i.e. VDR-BVs coincident with both #LCL eQTLs and GWAS SNPs [or in LD]). We're going to have a hard time over this with Ram, but the data justify this, IMHO.
#Re: LCL eQTLs
#Are there not tables of LCL eQTL SNPs cross-referenced with the genes whose expression they correlate with? If so, this is what I am interested in. I don't see why you #should recall eQTLs using LCL RNA-Seq data - surely this has already been done!?

#The sign denotes the direction of the nonreference allele (i.e. rvalue<0 means that nonreference allele has lower expression) 
#CONVERT THIS TO ANCESTRAL/DERIVED

#ultimately convert all these different geuvadis datasets into a unique bed with labels for ethnicity,gene,effect

#attention: indels have a decimal .5 in position 7, skip them

#typical input structure
#SNP_ID
#ID
#GENE_ID
#PROBE_ID
#CHR_SNP
#CHR_GENE
#SNPpos
#TSSpos
#distance
#rvalue
#pvalue
#log10pvalue

#rs3832000
#-
#ENSG00000142632.10
#ENSG00000142632.10
#1
#1
#16533832
#16539140
#5308
#-0.765767800472289
#3.68318824172032e-18
#17.433776084498

#1	SNP_ID : Variant identifier according to dbSNP137; position-based identifier for variants that are not in dbSNP (see Supplementary material pp 45) 
#2	ID : Null (-)
#3	GENE_ID : Gene identifier according to Gencode v12, miRBase v18, repeats based on their start site
#4	PROBE_ID : Quantitative trait identifier; the same as GENE_ID expect for:
#		Exons: GENEID_ExonStartPosition_ExonEndPosition
#		Transcript ratios: Transcript identifier according to Gencode v12
#5	CHR_SNP : Chromosome of the variant
#6	CHR_GENE : Chromosome of the quantitative trait
#7	SNPpos : Position of the variant
#8	TSSpos : Transcription start site of the gene/QT
#9	Distance : | SNPpos - TSSpos | 
#10	rvalue	: Spearman rank correlation rho (calculated from linear regression slope). The sign denotes the direction of the nonreference allele (i.e. rvalue<0 means that nonreference allele has lower expression) 
#11	pvalue : Association p-value
#12	log10pvalue : -log10 of pvalue

#bed output format 
#chr	pos-1	pos	rsid|geneid|rvalue|population|type
#1       1323143 1323144 rs141100746|ENSG00000242485.1|-0.265073362272798|CEU|exonQTL    2.05360647136187e-07

my $infile;
my $pop; #ceu/yri
my $type; #geneQTL/exonQTL/trratio
GetOptions(
        'i=s'      =>\$infile,
        'pop=s'    =>\$pop,
        'type=s'   =>\$type
);
if(!$infile){
	print "USAGE: do_funseq_GEUVADIS_to_bed.pl -i=<INFILE> \n";
    print "<INFILE> GEUVADIS GZIPPED tsv file\n";
    
    exit 1;
}
my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
my $output  = $directory . $basename . '.bed';

if($basename =~ /EUR/){
	$pop = 'CEU';
}elsif($basename =~ /YRI/){
	$pop = 'YRI';
}else{
	print STDERR "Error: impossible to determine ethnicity from filename: $basename. Aborting.\n";
	exit -1;
}
if($basename =~ /gene/){
	$type = 'geneQTL';
}elsif($basename =~ /exon/){
	$type = 'exonQTL';
}elsif($basename =~ /trratio/){
	$type = 'trratioQTL';
}else{
	print STDERR "Error: impossible to determine QTL type from filename: $basename. Aborting.\n";
	exit -1;
}

open (my $outstream,  q{>}, $output) or die("Unable to open $output : $!");
tie *FILE,   'IO::Zlib', $infile, "rb";
while (<FILE>)	{ 
	chomp;
	next if($_ =~ /^SNP_ID/);
	my @fields = split("\t", $_);
	
	#skip indels
	if($fields[6] =~ /\./){
		print STDERR "Warning: this appears to be an indel: $fields[6]. Skipping..\n";
		next;
	}
	#$5,$7-1 "",$7 "",$1"|"$3"|"$10"|CEU|trratioQTL",$11
	my $name_field = $fields[0] . '|' . $fields[2] . '|' . $fields[9] . '|' . $pop . '|' . $type;
	print $outstream $fields[4] . "\t" . ($fields[6]-1) . "\t" . $fields[6] . "\t" . $name_field . "\t" . $fields[10] . "\n"; 
}
close FILE;
close $outstream;

