#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#20/5/2014
#Takes the outfile from a SNPTEST analysis and produces a manhattan plot using R
#it uses qqman
#http://haldanessieve.org/2014/05/19/author-post-qqman-an-r-package-for-visualizing-gwas-results-using-q-q-and-manhattan-plots/
#http://cran.r-project.org/web/packages/qqman/

#first turns the input file in the requested input for qqman
#from:
#chr peak_id max_peak_rd alternate_ids rsid chromosome position alleleA alleleB index average_maximum_posterior_call info cohort_1_AA cohort_1_AB cohort_1_BB cohort_1_NULL all_AA all_AB all_BB all_NULL all_total all_maf missing_data_propo
#chr1 3 14.5761671611717 --- rs2488996 NA 974296 A G 1 1 1 2 6 16 0 2 6 16 0 24 0.208333 0 0.554959 0.413735 0.225809 NA


#to

#	SNP CHR BP         P
#1 rs1   1  1 0.9148060
#2 rs2   1  2 0.9370754
#3 rs3   1  3 0.2861395
#4 rs4   1  4 0.8304476
#5 rs5   1  5 0.6417455
#6 rs6   1  6 0.5190959

#note from qqman author:
#CHR column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again.


#todo remove cases with one GT < 1?

my $infile;
my $method;
my $ymax;
GetOptions(
        'i=s'         =>\$infile,
        'm=s'         =>\$method,
        'y=i'         =>\$ymax
);
if(!$infile){
     print "USAGE: do_manhattanplotR_from_association.pl -i=<ASSOCIATIONS> -m=<METHOD> -y=<YMAX>\n";
     print "<ASSOCIATIONS> file containing global results from SNPTEST\n";
     print "<METHOD> one of [pval|qval|bf]\n";
     print "<YMAX> if present, overrides the automatic setting to the highest BF\n";
     exit 1;
}
if(!$method){
     print "USAGE: do_manhattanplotR_from_association.pl -i=<ASSOCIATIONS> -m=<METHOD> -y=<YMAX>\n";
     print "<ASSOCIATIONS> file containing global results from SNPTEST\n";
     print "<METHOD> one of [pval|qval|bf]\n";
     print "<YMAX> if present, overrides the automatic setting to the highest BF\n";
     exit 1;
}
unless(  ($method eq 'qval')  or ($method eq 'pval') or ($method eq 'bf')  ){
     print "USAGE: do_manhattanplotR_from_association.pl -i=<ASSOCIATIONS> -m=<METHOD> -y=<YMAX>\n";
     print "<ASSOCIATIONS> file containing global results from SNPTEST\n";
     print "<METHOD> one of [pval|qval|bf]\n";
     print "<YMAX> if present, overrides the automatic setting to the highest BF\n";
     exit 1;
}

my $PATH_CODE = '/net/isi-scratch/giuseppe/tools/R-3.0.1/bin/Rscript';
my($basename, $directory) = fileparse($infile);
$basename =~ s/(.*)\..*/$1/;
my $temp_file   = $directory   . 'manhattan_' . $basename . '.tmp';
my $temp_file_R = $directory   . 'manhattan_' . $basename . '.R';
my $temp_file_PDF = $directory . 'manhattan_' . $basename . '.pdf';
my $temp_file_SVG = $directory . 'manhattan_' . $basename . '.svg';


#build qqman input file
open (my $instream,     q{<}, $infile) or die("Unable to open $infile : $!");
my %plot_structure;
while(<$instream>){	
	chomp;
	next if($_ =~ /comment/);# header
	my ($chr, $rs_id, $bp, $value );
	
	if($method eq 'pval'){
		($chr, $rs_id, $bp, $value ) = (split ' ')[0,4,6,23];
	}elsif($method eq 'qval'){
		($chr, $rs_id, $bp, $value ) = (split ' ')[0,4,6,28];
	}elsif($method eq 'bf'){
		($chr, $rs_id, $bp, $value ) = (split ' ')[0,4,6,23]; #actually same as pval
	}
	#convert chr string
	if($chr =~ /chr(.+)/){
		$chr = $1;
	}
	$plot_structure{$chr}{$bp}{$rs_id} = $value;
}
close $instream;

open (my $tempstream,     q{>}, $temp_file) or die("Unable to open $temp_file : $!");

if($method eq 'bf'){
	print $tempstream "SNP\tCHR\tBP\tBF\n";
}else{
	print $tempstream "SNP\tCHR\tBP\tP\n";
}
foreach my $chr (sort {$a<=>$b} keys %plot_structure){
	foreach my $bp ( sort {$a<=>$b} keys %{ $plot_structure{$chr} } ){
		foreach my $rs_id ( sort  keys %{ $plot_structure{$chr}{$bp} } ){
			print $tempstream $rs_id . "\t" . $chr . "\t" . $bp . "\t" .  $plot_structure{$chr}{$bp}{$rs_id} . "\n";
		}
	}
}
close $tempstream;

#create R script and call R
open ($tempstream,     q{>}, $temp_file_R) or die("Unable to open $temp_file_R : $!");

if($method eq 'bf'){
	print $tempstream 'source("/net/isi-scratch/giuseppe/tools/qqman/qqman_bayesfactors/R/manhattan.R")' . "\n";
}else{
	print $tempstream 'library(qqman)' . "\n";
}
print $tempstream 'data<-read.delim(\'' . $temp_file . '\',header=T)' . "\n";
#print $tempstream 'pdf(\'' . $temp_file_PDF . '\')' . "\n";
print $tempstream 'svg(\'' . $temp_file_SVG . '\')' . "\n";
#eg
#manhattan(gwasResults, col = c("blue4", "orange3"), main = "Results from simulated trait", genomewideline = FALSE, suggestiveline = FALSE)
#print $tempstream 'manhattan(subset(data, CHR == 6))' . "\n";
#cex = 0.5
#ymax = 8

if($ymax){
	if ($method eq 'pval'){
		print $tempstream 'manhattan(data, col = c("black","#666666","#CC6600"), ymax=' . $ymax  . ', main = "Manhattan plot of VDR-QTL p-values",  genomewideline = FALSE, suggestiveline = FALSE)' . "\n";
		print $tempstream 'qq(data$P, main = "Q-Q plot of VDR-QTL p-values")' . "\n";
	}elsif($method eq 'qval'){
		print $tempstream 'manhattan(data, col = c("black","#666666","#CC6600"), ymax=' . $ymax  . ', main "Manhattan plot of VDR-QTL q-values",  genomewideline = FALSE, suggestiveline = FALSE)' . "\n";
		print $tempstream 'qq(data$P, main = "Q-Q plot of VDR-QTL p-values")' . "\n";
	}elsif($method eq 'bf'){
		print $tempstream 'manhattan(data, col = c("black","#666666","#CC6600"), ymax=' . $ymax  . ', main = "Manhattan plot of VDR-QTL Bayes Factors (SNPTEST)",  genomewideline = FALSE, suggestiveline = FALSE)' . "\n";
	}
}else{
	if ($method eq 'pval'){
		print $tempstream 'manhattan(data, col = c("black","#666666","#CC6600"), main = "Manhattan plot of VDR-QTL p-values",  genomewideline = FALSE, suggestiveline = FALSE)' . "\n";
		print $tempstream 'qq(data$P, main = "Q-Q plot of VDR-QTL p-values")' . "\n";
	}elsif($method eq 'qval'){
		print $tempstream 'manhattan(data, col = c("black","#666666","#CC6600"), main = "Manhattan plot of VDR-QTL q-values",  genomewideline = FALSE, suggestiveline = FALSE)' . "\n";
		print $tempstream 'qq(data$P, main = "Q-Q plot of VDR-QTL p-values")' . "\n";
	}elsif($method eq 'bf'){
		print $tempstream 'manhattan(data, col = c("black","#666666","#CC6600"), main = "Manhattan plot of VDR-QTL Bayes Factors (SNPTEST)",  genomewideline = FALSE, suggestiveline = FALSE)' . "\n";
	}	
}
print $tempstream 'dev.off()' . "\n";
close $tempstream;


system "$PATH_CODE $temp_file_R";
unlink $temp_file;
unlink $temp_file_R;
