#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

#27/6/2016
#script to convert Uniprot IDs to Ensembl Protein IDs in Rebecca's inputs
#This is the second step, and must be used after the first script that cleans up
#the Uniprot Ids then those IDs need to be queried online (I used the Uniprot ID converter)



#input
#1) 

#csv table from rebecca as follows

#Protein	Description	MSS11223 SINQ Co culture MIC SIn	MSS11224 SINQ Myocytes MIC SIn	ratio Co/Myo	fold regulation Co/myo
#sp|Q8VHU4|ELP1_RAT	Elongator complex protein 1 OS=Rattus norvegicus GN=Ikbkap PE=2 SV=1		3.39E-008	0	#DIV/0!

#2) mapping table uniprot > ensembl obtained eg from Uniprot
#From    To
#A0A096MIU5      ENSRNOP00000054552


#clean id from field 0
#query it against hash uniprot > ensembl
#save field 0, new id, field 2,3,4
#for now, skip all those that have DIV in field 4
#for the future,
#if 2 or 3 are empty then use the checks Rebecca wants


#output
#tsv table with
#ORIG_ID\tENSP_ID\tGENE_NAME\tRATIO


#This version of the script uses the Ensemb API to get the gene name

my $infile;
my $map;
GetOptions(
    'i=s'  =>\$infile,
    'm=s'  =>\$map
);


my $USAGE = "USAGE: $0 -i=<INFILE> -m=<MAPFILE>\n" . 
            "<INFILE> csv obtained from Rebecca's tables including expression level column\n" .
            "<MAPFILE> uniprot to ensembl prot map obtained from uniprot webservice\n";
die $USAGE if(!$infile || !$map);

my $registry = 'Bio::EnsEMBL::Registry';
$registry -> load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

#my $gene_adaptor = $registry->get_adaptor('Rattus Norvegicus', 'core', 'Gene');
my $translation_adaptor = $registry->get_adaptor( "Rattus Norvegicus", "core","translation" );
die "Unable to obtain translation adaptor" unless($translation_adaptor);


#----------
#1. Get map
#----------
my %uniprot_to_ensembl;
open (my $instream,  q{<}, $map) or die("Unable to open $map : $!");
while(<$instream>){
    chomp;
    my ($id_uniprot, $id_ens) = split("\t", $_);
    $uniprot_to_ensembl{$id_uniprot} = $id_ens;
}
close $instream;


my $count=0;
print STDOUT "ORIG_ID\tENSP_ID\tENSG_NAME\tEXPR_RATIO\n";
open ($instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^Protein/); #header
	my $translation; my $gene_name;
	my $index;
    my $INPUT_ID; my $ENS_ID;
	my ($input_id_string, $treat, $control, $ratio_tr_ctrl) = (split /\t/)[0,2,3,4];

    if($input_id_string =~ /.*\|(.*)\|.*/){
        $INPUT_ID = $1;
	}else{
        print STDERR "Warning: unable to extract UniprotID from this $input_id_string. Skipping..\n";
        next;
    }

    #remove any isoforms 
    $INPUT_ID =~ s/(.*)\-\d+/$1/;
    $ENS_ID = $uniprot_to_ensembl{$INPUT_ID};

    if(!$ENS_ID){
        print STDERR "No ID mapping available for $INPUT_ID. Skipping..\n";
        next
    }

    #get gene name from translation id
    $translation = $translation_adaptor->fetch_by_stable_id($ENS_ID);
    $gene_name = $translation->transcript->get_Gene->external_name();
	$gene_name = '' if(!$gene_name);

	print STDOUT $INPUT_ID  . "\t" . $ENS_ID . "\t" . $gene_name . "\t" . $ratio_tr_ctrl . "\n" if($ENS_ID);
}
close $instream;
print STDERR 'count missed: ' . $count . "\n";
