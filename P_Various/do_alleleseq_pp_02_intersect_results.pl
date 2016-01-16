#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;


#I want to obtain statistics for the intersection of the allele specific events 
#I will provide this script with a directory containing all interestingHets.txt files

#For each file, for each line get a unique key for the line (chr-pos) and put it in a hash with a counter
#Note that a chr-pos can map to more than one rsID. Take all of them?

#output
#a CONSENSUS interestinghets.txt file with all chr-pos events appearing in at least 2 samples. 
#a chr-pos -> number_of_observations to build a histogram of counts

#The objective would be to build a ggplot2 matrix plot, like a heatmap, indicating which events are there.

my $dir;
my $allsites; #flag. If set, gets ASB sites even when they don't fall in a peak bed interval
my $MIN_OVERLAP; #minimum overlap for printing.
GetOptions(
        'd=s'  => \$dir,
        'allsites'   =>\$allsites,
        'o=i' =>\$MIN_OVERLAP
);

if(!$dir){
     print "USAGE: do_alleleseq_intersect_results.pl -d=<DIR> -o=<MIN_SAMPLE_OVERLAP> -allsites\n";
     print "<DIR> directory containing the input files in the format <interestingHets_NAxxxxx.txt>\n";
     print "(opt)<MIN_SAMPLE_OVERLAP> minimum overlap for printing. eg -o=2 means only asb events overlapping at least 2 samples will be considered. Default:all\n";
     print "(opt)<allsites> flag; if set, retrieve asb sites even when they did not fall in peak interval file.\n";
     exit 1;
}

$MIN_OVERLAP = 1 if(!$MIN_OVERLAP);
if($MIN_OVERLAP < 1){
	print "ERROR: $MIN_OVERLAP is not a positive integer\n";
	exit -1;
}

chdir $dir;
my @files = <interestingHets_NA*.txt>;

my %asymevent_to_sample;
my %asymevent_to_count;
my %asymevent_is_peak;
my %sample_list;
my %hash_fullentry;
my $instream;
my $header;

my $output_consensus;
my $output_matrix;
if($allsites){
	$output_consensus = "interestingHets_consensus_allsites_o" .  $MIN_OVERLAP . '.txt';
	$output_matrix = "interestingHets_matrix_allsites_o" .  $MIN_OVERLAP . '.txt';
}else{
	$output_consensus = "interestingHets_consensus_inpeaks_o" .  $MIN_OVERLAP . '.txt';
	$output_matrix = "interestingHets_matrix_inpeaks_o" .  $MIN_OVERLAP . '.txt';	
}


#chrm	snppos   	ref	mat_gtyp	pat_gtyp	c_gtyp	phase	mat_all	pat_all	cA	cC	cG	cT	winning	SymCls	SymPval	BindingSite	cnv
#1	1080920	G	S	W	R	PHASED	G	A	2	0	7	0	M	Sym	0.1796875	1	1.0

foreach my $file (@files){
	my $sample_id = $file;
	$sample_id =~ s/interestingHets_(NA\d{5}).*/$1/;
	$sample_list{$sample_id} = 1;
	open ($instream,  q{<}, $file) or die("Unable to open $file : $!");
	while(<$instream>){
		chomp;
		next if($_ eq '');
		if($_ =~ /^chrm/){
			$header = $_;
			next;
		}
		my ($a_chr,$a_snppos,$a_ref,$a_symCls,$a_symPval,$a_bindingSite) = (split /\t/)[0,1,2,14,15,16];
		next if(!$a_chr);
		next if(!$a_snppos);
		next if(!$a_ref);
		next if(!$a_symCls);
		next if(!$a_symPval);
		next unless ($a_symCls =~ /Asym/);
		
		if(!$allsites){
			next unless ($a_bindingSite =~ /1/);
		}
		
		my $event_key  = $a_chr . '-' . $a_snppos . '-' . $a_ref;
		my $full_entry = $a_chr . '-' . $a_snppos . '-' . $a_ref . '-' . $a_bindingSite;
		
		
		$asymevent_is_peak{$event_key} = $a_bindingSite;
		$asymevent_to_sample{$event_key}{$sample_id} = 1;
		if(!$asymevent_to_count{$event_key}){
			$asymevent_to_count{$event_key} = 1;
		}else{
			$asymevent_to_count{$event_key} += 1;
		}
	}
	close $instream;	
}

#output1: consensus interestingHet.txt file, with additional column indicating how many times the event is observed across the samples  (including 1 time) and which samples have it
#I want to output the following columns

open (my $outstream,  q{>}, $output_consensus) or die("Unable to open $output_consensus : $!");
print $outstream "chrm\tsnppos\tref\tBindingSite\toverlap\tsampleIDS" . "\n";
foreach my $event (sort keys %asymevent_to_sample){
	if($asymevent_to_count{$event} > ($MIN_OVERLAP - 1) ){ 
		my @samples_for_this_event;
		my ($chr,$snppos,$a_ref) = split('-', $event);
		foreach my $sample (sort keys %{ $asymevent_to_sample{$event} }){
			push(@samples_for_this_event, $sample);
		}
		my $line = $chr . "\t" . 
					$snppos . "\t" . 
					$a_ref . "\t" . 
					$asymevent_is_peak{$event} . "\t" . 
					$asymevent_to_count{$event} . "\t" . 
					join(',', @samples_for_this_event) . "\n";
		print $outstream $line;
	}
} 
close $outstream;



#output2: matrix: asb event X sample, for events seen in at least 2 samples

open ($outstream,  q{>}, $output_matrix) or die("Unable to open $output_matrix : $!");

#header: event,sample names
my @all_samples;
foreach my $sample (sort keys %sample_list){ push(@all_samples, $sample); }
print $outstream 'ASB_event' . "\t" . join("\t", @all_samples) . "\n";

foreach my $event (sort keys %asymevent_to_sample){
	
	if($asymevent_to_count{$event} > ($MIN_OVERLAP - 1) ){
		#need to create a string of 1s and 0s based on the presence or absence of a sample in the full list
		my @samples_for_this_event;
		foreach my $sample (sort keys %sample_list){
			if($asymevent_to_sample{$event}{$sample}){
				push(@samples_for_this_event, '1');
			}else{
				push(@samples_for_this_event, '0');
			}	
		}
		 

#		foreach my $sample (sort keys %{ $asymevent_to_sample{$event} }){
#			if($asymevent_to_sample{$event}{$sample}){
#				push(@samples_for_this_event, '1');
#			}else{
#				push(@samples_for_this_event, '0');
#			}	
#		}
		
		print $outstream $event . "\t" . join("\t", @samples_for_this_event) . "\n";
	}
}
close $output_matrix;