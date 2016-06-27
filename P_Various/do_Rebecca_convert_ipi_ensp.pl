#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

#old script to convert Rebecca's IPI ids to ENSP ids using the Ensembl API

#input
#csv table from rebecca as follows
#Accession #	Expression
#IPI00231783.5	0.001
#IPI00776952.1	0.002
#IPI00231699.5	0.003
#IPI00211053.6	0.003
#IPI00780332.1	0.006
#IPI00393489.2	0.006
#IPI00189811.1	0.007
#IPI00202616.1	0.008
#IPI00212908.2	0.008
#IPI00213057.3	0.008
#IPI00198497.1	0.009
#IPI00951953.1	0.009

#output
#tsv table with
#IPI\tENSP_ID\tSHR_WKY_expression_value

#Idea: query ensembl with each ipi,GET the ensp for each line
#NOTA ad ogni ipi puo corrispondere piu di un ensembl id

my $infile;
GetOptions(
        'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_data_Rebecca.pl -i=<INFILE>\n";
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



