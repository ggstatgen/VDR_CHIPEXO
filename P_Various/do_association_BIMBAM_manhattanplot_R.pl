#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#run this after do_association_BIMBAM.pl to get a manhattan plot of the results

#26/5/2014
#Takes the outfile from a BIMBAM analysis and produces a manhattan plot using R
#it uses qqman
#http://haldanessieve.org/2014/05/19/author-post-qqman-an-r-package-for-visualizing-gwas-results-using-q-q-and-manhattan-plots/
#http://cran.r-project.org/web/packages/qqman/

#first turns the input file in the requested input for qqman
#from:
#chr     pos     rsnum   bf      pval    mu      a       d
#1       713979  rs117217250     +0.000   4.88e-01       1.687   0.000   0.000

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

my $infile;
my $ymax;
GetOptions(
        'i=s'         =>\$infile,
        'y=i'         =>\$ymax
);
if(!$infile){
     print "USAGE: do_association_BIMBAM_manhattanplot_R.pl -i=<ASSOCIATIONS> -y=<YMAX>\n";
     print "<ASSOCIATIONS> file containing global results from BIMBAM\n";
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
my $temp_file_PNG = $directory . 'manhattan_' . $basename . '.png';

#build qqman input file
open (my $instream,     q{<}, $infile) or die("Unable to open $infile : $!");
my %plot_structure;
while(<$instream>){	
	chomp;
	next if($_ =~ /^chr/);# header
	my ($chr, $bp, $rs_id, $value ) = (split "\t")[0,1,2,3];
	
	#convert chr string
	if($chr =~ /chr(.+)/){
		$chr = $1;
	}
	$plot_structure{$chr}{$bp}{$rs_id} = $value;
}
close $instream;

open (my $tempstream,     q{>}, $temp_file) or die("Unable to open $temp_file : $!");
print $tempstream "SNP\tCHR\tBP\tBF\n";

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
print $tempstream 'source("/net/isi-scratch/giuseppe/tools/qqman/qqman_bayesfactors/R/manhattan.R")' . "\n";
print $tempstream 'data<-read.delim(\'' . $temp_file . '\',header=T)' . "\n";
#print $tempstream 'pdf(\'' . $temp_file_PDF . '\')' . "\n";
#print $tempstream 'png(\'' . $temp_file_PNG . '\')' . "\n";
print $tempstream 'svg(\'' . $temp_file_SVG . '\')' . "\n";
#eg
#manhattan(gwasResults, col = c("blue4", "orange3"), main = "Results from simulated trait", genomewideline = FALSE, suggestiveline = FALSE)
#print $tempstream 'manhattan(subset(data, CHR == 6))' . "\n";
#cex = 0.5
#ymax = 8

if($ymax){
	#print $tempstream 'manhattan(data,cex=0.5,ymax=' . $ymax  . ', main = "Manhattan plot of VDR-QTL Bayes Factors (BIMBAM)",  genomewideline = FALSE, suggestiveline = FALSE)' . "\n";
	print $tempstream 'manhattan(data, col = c("black","#666666","#CC6600"), ymax=' . $ymax  . ', main = "Manhattan plot of VDR-QTL Bayes Factors (BIMBAM)",  genomewideline = FALSE, suggestiveline = FALSE)' . "\n";
	#print $tempstream 'manhattan(data, col = "#666666", ymax=' . $ymax  . ', main = "Manhattan plot of VDR-QTL Bayes Factors (BIMBAM)",  genomewideline = FALSE, suggestiveline = FALSE)' . "\n";
}else{
	#print $tempstream 'manhattan(data, cex=0.5, main = "Manhattan plot of VDR-QTL Bayes Factors (BIMBAM)",  genomewideline = FALSE, suggestiveline = FALSE)' . "\n";
	print $tempstream 'manhattan(data, col = c("black","#666666","#CC6600"), main = "Manhattan plot of VDR-QTL Bayes Factors (BIMBAM)",  genomewideline = FALSE, suggestiveline = FALSE)' . "\n";
	#print $tempstream 'manhattan(data, col = "#777777", main = "Manhattan plot of VDR-QTL Bayes Factors (BIMBAM)",  genomewideline = FALSE, suggestiveline = FALSE)' . "\n";
}
print $tempstream 'dev.off()' . "\n";
close $tempstream;

system "$PATH_CODE $temp_file_R";
unlink $temp_file;
unlink $temp_file_R;
