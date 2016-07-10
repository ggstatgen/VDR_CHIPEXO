#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max);


#TODO get rid of all the DIV rows for now

#27/6/2016
#This script DOES NOT actually convert the UNIPROT KB IDs
#Maybe String and Expressence can deal with them
#So this just cleans up the file, gets the gene name and the fold ration

#input
#csv table from rebecca as follows
#Protein	Description	MSS11223 SINQ Co culture MIC SIn	MSS11224 SINQ Myocytes MIC SIn	ratio Co/Myo	fold regulation Co/myo
#sp|Q8VHU4|ELP1_RAT	Elongator complex protein 1 OS=Rattus norvegicus GN=Ikbkap PE=2 SV=1		3.39E-008	0	#DIV/0!



#output
#tsv table with
#ORIG_ID\tGENE_NAME\tRATIO

my $infile;
GetOptions(
    'i=s'  =>\$infile
);


my $USAGE = "USAGE: $0 -i=<INFILE>\n" . 
            "<INFILE> csv obtained from Rebecca's tables including expression level column\n";
die $USAGE if(!$infile);


my $count=0;
#print STDOUT "ORIG_ID\tNAME\tEXPR_RATIO\n";
print STDOUT "ORIG_ID\tNAME\tTREAT\tCTRL\n";
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^Protein/); #header
    next if($_ =~ /\#DIV/); #either treatment or control are empty

    my $gene_name;
	my $index;
    my $INPUT_ID;
	my ($input_id_string, $treat, $control, $ratio_tr_ctrl) = (split /\t/)[0,2,3,4];

    if($input_id_string =~ /.*\|(.*)\|(.*)\_RAT/){
        $INPUT_ID = $1;
        $gene_name = $2;
	}else{
        print STDERR "Warning: unable to extract UniprotID and gene name from this $input_id_string. Skipping..\n";
        next;
    }
    #remove any isoforms 
    $INPUT_ID =~ s/(.*)\-\d+/$1/;
    
    #print STDOUT $INPUT_ID  . "\t" . ucfirst($gene_name) . "\t" . $ratio_tr_ctrl . "\n";
    print STDOUT $INPUT_ID  . "\t" . ucfirst($gene_name) . "\t" . $treat . "\t" . $control . "\n";

}
close $instream;
print STDERR 'count missed: ' . $count . "\n";
