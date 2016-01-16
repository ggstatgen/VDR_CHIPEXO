#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX 'log10';

#10/2/2015
#I found this way to do hinton plots
#http://tex.stackexchange.com/questions/155291/generate-a-hinton-diagram-using-pgfplots
#It uses a set of latex packages and I was able to make the code work on windows (you will need to run the code on the windows machine)

#This script simply converts the matrix form data produced by 
#the R PMCMR and used to test the phastcons scores:
#/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/PLOT_VDRRXR_phastcons
#into a linear format of the kind:
#0 0 2.766580
#1 0 -0.079900

#you could produce -log_10(pval) or maybe plot the chi-square. Maybe both

#INPUT: matrix obtained from R, in the following form:
#"","p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","p13","p14"
#"p2",0.999999929835591,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA
#...
#"p15" ...

#output
#v_p1	v_p2	[-log_10(p)|chisquare]

#ATM I replace NA with zeros, though I could probably also try to make the matrix symmetrical
#I add 1	1	0 to the first line to make the number or rows in the final plot right


#COPY THIS OUTPUT TO THE WINDOWS MACHINE and GENERATE PLOT WITHIN LATEX
my $infile;
my $type;
GetOptions(
        'i=s'	=>\$infile,
        't=s'	=>\$type
);
#$infile = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/PLOT_VDRRXR_phastcons/kruskal_matrix_pa_stats.txt";
#$type='stats';

my $SIGNIFICANCE_LEVEL = 0.05;

if(!$infile){
	print "USAGE: do_hinton_latexdata_from_kruskal -i=<INFILE> -t=<pvals|stats>\n";
	print "<INFILE> matrix of pvals or stats produced using R PMCMR package\n";
	print "<TYPE> one of [pvals|stats], whether <INFILE> contains pvals or test statistics";
	print "NOTE1: when type=pvals, values are log10 transformed before producing output\n";
	print "NOTE2: significance level for plotting TYPE=pvals set at $SIGNIFICANCE_LEVEL. Change in script if desired.\n";
    exit 1;
}
if(!$type){
	print "USAGE: do_hinton_latexdata_from_kruskal -i=<INFILE> -t=<pvals|stats>\n";
	print "<INFILE> matrix of pvals or stats produced using R PMCMR package\n";
	print "<TYPE> one of [pvals|stats], whether <INFILE> contains pvals or test statistics";
	print "NOTE1: when type=pvals, values are log10 transformed before producing output\n";
	print "NOTE2: significance level for plotting TYPE=pvals set at $SIGNIFICANCE_LEVEL. Change in script if desired.\n";
    exit 1;
}

print STDOUT "1\t1\t0\n";

my @X_array;
my %positional_hash;
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	my $index_X; my $index_Y;
	my @data = split(',',$_);
	my $matrix_elements = scalar(@data); #the number of positions in the motif
	my $element = shift @data;
	
	next if($element eq "\"\"");

	if($element =~ /\"p(\d+)\"/){
		$index_Y = $1;
	}else{
		print "Y coordinate not recognised: $element. Aborting..\n";
		exit -1;
	}
	
	$index_X = 1;
	if($type eq 'stats'){#----------------------------------------------------
		foreach my $value (@data){
			$value = 0 if($value eq 'NA');	
			print STDOUT $index_X, "\t", $index_Y, "\t", $value, "\n"; 	
			$index_X++;
			if($index_X == 15){
				print STDOUT $index_X, "\t", $index_Y, "\t", '0', "\n";
			}	
		}		
	}elsif($type eq 'pvals'){#------------------------------------------------
		foreach my $value (@data){
			if($value eq 'NA'){
				$value = 0;
			}else{
				my $SIGNIFICANT = 1 if($value <= $SIGNIFICANCE_LEVEL);
				#trick to get the latex program on the laptop to plot the significant in different colour than non significant.
				#latex program uses sign to distinguish colour so I will artificially add sign '-' to all signficant
			
				$value = -log10($value);
				$value = 0 if($value eq '-0');	
				$value = '-' . $value if($SIGNIFICANT);				
			}
			print STDOUT $index_X, "\t", $index_Y, "\t", $value, "\n"; 	
			$index_X++;
			if($index_X == 15){
				print STDOUT $index_X, "\t", $index_Y, "\t", '0', "\n";
			}	
		}
	}else{#------------------------------------------------------------------
		print STDERR "ERROR: type: $type is not recognised. Aborting..\n";
		exit -1;
	}
}
close $instream;
