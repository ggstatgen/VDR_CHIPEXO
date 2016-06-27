#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

#27/6/2016
#script to clean up data file from Holger/Rebecca of anything to necessary and at the
#the same time convert Uniprot to EnsP using the Ensembl API

#The idea is that you first run this, then for all the unsolved cases, you try biomart.


#input
#csv table from rebecca as follows

#Protein	Description	MSS11223 SINQ Co culture MIC SIn	MSS11224 SINQ Myocytes MIC SIn	ratio Co/Myo	fold regulation Co/myo
#sp|Q8VHU4|ELP1_RAT	Elongator complex protein 1 OS=Rattus norvegicus GN=Ikbkap PE=2 SV=1		3.39E-008	0	#DIV/0!

#clean id from field 0
#query it against Ensembl API
#save field 0, new id, field 2,3,4
#for now, skip all those that have DIV in field 4
#for the future,
#if 2 or 3 are empty then use the checks Rebecca wants


#output
#tsv table with
#ORIG_ID\tENSP_ID\tRATIO


my $infile;
GetOptions(
        'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: $0 -i=<INFILE>\n";
     print "<INFILE> csv obtained from Rebecca's tables including expression level column\n";
     exit 1;
}

my $registry = 'Bio::EnsEMBL::Registry';
$registry -> load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

#my $gene_adaptor = $registry->get_adaptor('Rattus Norvegicus', 'core', 'Gene');
my $translation_adaptor = $registry->get_adaptor( "Rattus Norvegicus", "core","translation" );
die "Unable to obtain translation adaptor" unless($translation_adaptor);

my $count=0;
print "ORIG_ID\tENSP_ID\tENSG_NAME\tEXPR\n";
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /Accession/); #header
	my $gene_name;
	my $index;
	my ($INPUT_ID, $expression_value) = (split /\t/)[0,1];
	#print $INPUT_ID, "\t", $expression_value, "\n";
	my $ENS_ID;
	
	#refseq translations
	my @translations = @{$translation_adaptor->fetch_all_by_external_name($INPUT_ID)};
	if(!@translations){
		#print STDERR "\nNo translations for $INPUT_ID\n";
		$count++;
		next;
	}
		
	#now we suppose there are >1 translations. Let's get the most recent gene name:
#	my @genes = @{$gene_adaptor->fetch_all_by_external_name($INPUT_ID)};
#	next if(!@genes);#should never get in here - if it does the translation is some random crap
#	if (scalar @genes > 1){ #I'll choose the last version
#		my $max;my $x = 0;my @versions;$index = 0;
#		foreach my $gene(@genes) { push(@versions, $gene->version()); } 
#		for ( @versions ){
#			if (!defined $max or $_ > $max){ $index = $x; $max = $_; }
#    		$x++;
#		}
#		$gene_name = $genes[$index]->external_name();
#	}else{#only one gene object for this IPI
#		$gene_name = $genes[0]->external_name();
#	}
	
	#same for translations now
	if(scalar @translations > 1){ #I'll chose the last version
		my $max;my $x = 0;my @versions;$index = 0;
		foreach my $translation (@translations) { push(@versions, $translation->version()); }
		for ( @versions ){
			if (!defined $max or $_ > $max){ $index = $x; $max = $_; }
    			$x++;
		}
		$ENS_ID = $translations[$index]->stable_id();
		#get the gene name # SUGGESTION DAN STAINES
		my $transcript = $translations[$index]->transcript();
		my $gene = $transcript->get_Gene();
		$gene_name = $gene->external_name();
		if(!$gene_name){
			$gene_name = ($gene->display_xref()) ? $gene->display_xref()->primary_id() : $gene->stable_id();
		}
	}else{#only one translation object for this INPUT ID
		$ENS_ID = $translations[0]->stable_id();
		my $transcript = $translations[0]->transcript();
		my $gene = $transcript->get_Gene();
		$gene_name = $gene->external_name();
		if(!$gene_name){
			$gene_name = ($gene->display_xref()) ? $gene->display_xref()->primary_id() : $gene->stable_id();
		}
	}
	$gene_name = '' if(!$gene_name);
	print $INPUT_ID  . "\t" . $ENS_ID . "\t" . $gene_name . "\t" .  $expression_value . "\n" if($ENS_ID);
}
close $instream;
#print STDERR 'count: ' . $count . "\n";
	

	
	#print $INPUT_ID . "\n";
#	#my $gene = $gene_adaptor->fetch_by_display_label($INPUT_ID);
#	my @genes = @{$gene_adaptor->fetch_all_by_external_name($INPUT_ID)};
#	foreach my $gene (@genes){
#		
#	print "GENE ", $gene->stable_id(), "\n";
#	print_DBEntries( $gene->get_all_DBEntries() );
#	
#	foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
#    		print "TRANSCRIPT ", $transcript->stable_id(), "\n";
#    		print_DBEntries( $transcript->get_all_DBEntries() );
#
#    		# Watch out: pseudogenes have no translation
#    		if ( defined $transcript->translation() ) {
#        		my $translation = $transcript->translation();
#        		print "TRANSLATION ", $translation->stable_id(), "\n";
#        		print_DBEntries( $translation->get_all_DBEntries() );
#    		}
#	}
#	
#	}


#sub print_DBEntries{
#    my $db_entries = shift;

    #foreach my $dbe ( @{$db_entries} ) {
    #    printf "\tXREF %s (%s)\n", $dbe->display_id(), $dbe->dbname();
    #}
#}



