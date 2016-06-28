#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

#27/6/2016
#script to clean up data file from Holger/Rebecca of anything to necessary and return a list of 
#UniprotIDs to send to the following interface:
#http://www.uniprot.org/uploadlists/
#for conversion to ENSP

#input
#csv table from rebecca as follows

#Protein	Description	MSS11223 SINQ Co culture MIC SIn	MSS11224 SINQ Myocytes MIC SIn	ratio Co/Myo	fold regulation Co/myo
#sp|Q8VHU4|ELP1_RAT	Elongator complex protein 1 OS=Rattus norvegicus GN=Ikbkap PE=2 SV=1		3.39E-008	0	#DIV/0!

#clean id from field 0 and return it



my $infile;
GetOptions(
    'i=s'  =>\$infile,
);
if(!$infile){
    print "USAGE: $0 -i=<INFILE>\n";
    print "<INFILE> csv obtained from Rebecca's tables\n";
    exit 1;
}

open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^Protein/); #header
	my $gene_name;
	my $index;
	my ($input_id_string) = (split /\t/)[0];
    my $INPUT_ID; 

    if($input_id_string =~ /.*\|(.*)\|.*/){
        $INPUT_ID = $1;
	}else{
        print STDERR "Warning: unable to extract UniprotID from this $input_id_string. Skipping..\n";
        next;
    }
    $INPUT_ID =~ s/(.*)\-\d+/$1/;
    print $INPUT_ID, "\n"; next;
}
close $instream;
