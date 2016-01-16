#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#I have a set of files from my htseq count, one per sample, one row per peak. I want a matrix as follows:
#      sample1	sample2	...	sampleN
#peak1	...	...	...
#peak2	...	...	...
#...
#peakN

#I will then pass this DEseq2 or Limma for normalisation
#INPUT: directory with all the count files to process
#output text file with matrix of counts

my $in_counts_dir;
my $resultdir = 'd_out';
my $PEAKS;
GetOptions(
        'counts=s'       =>\$in_counts_dir#,
        #'peaknumber=i'   =>\$PEAKS
);
if(!$in_counts_dir){
	print "USAGE: do_htseq_get_matrix_from_counts.pl -counts=<HTSEQ_COUNT_DIR>\n";
	print "<HTSEQ_COUNT_DIR> directory containing all the count text files for the required samples\n";
	print "WARNING: this script assumes the peak id returned by the peak caller IS A NUMBER (eg GEM)\n";
    exit 1;
}
#if(!$in_counts_dir){
#	print "USAGE: do_htseq_get_matrix_from_counts.pl -counts=<HTSEQ_COUNT_DIR> -peaks=<NUM_PEAKS>\n";
#	print "<HTSEQ_COUNT_DIR> directory containing all the count text files for the required samples\n";
#	print "<PEAKS> total number of peaks used for the Htseq count\n";
#    exit 1;
#}
#if(!$PEAKS){
#	print "USAGE: do_htseq_get_matrix_from_counts.pl -counts=<HTSEQ_COUNT_DIR> -peaks=<NUM_PEAKS>\n";
#	print "<HTSEQ_COUNT_DIR> directory containing all the count text files for the required samples\n";
#	print "<PEAKS> total number of peaks used for the Htseq count\n";
#    exit 1;
#}
chdir $in_counts_dir;
my @files = <*_counts.txt>;

my %count_matrix;
my %sample_dictionary;
my %peak_dictionary;
foreach my $infile (@files) {
	my $id;
	#current format VDR_GM12264_reads_htseq_counts.txt
	if($infile =~ /\_GM(\d+)\_/){
		 $id = 'NA' . $1;
		 $sample_dictionary{$id} = 1;
	}
	#preinitialise matrix with zeros
	#foreach my $i ( 0..($PEAKS-1) ){ $count_matrix{$i}{$id} = 0; }
		
	open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
	while(<$instream>){
		chomp;
		next if($_ =~ /^\_/); # _mapped, _ambiguous, _no_feature, etc
		my ($peak_ID, $count) = (split / /)[0,1];
		$peak_dictionary{$peak_ID} = 1;
		$count_matrix{$peak_ID}{$id} = $count;
	}
	close $instream;
}

#cycle through all peak ids and all samples and fill the matrix with zeros when a count is missing
foreach my $sample (sort keys %sample_dictionary){
	foreach my $peak (sort keys %peak_dictionary){
		$count_matrix{$peak}{$sample} = 0 if(!$count_matrix{$peak}{$sample});
	}
}

#print header
#peakID sample1 sample2 ... sampleN
print 'peakID';
foreach my $sample (sort keys %sample_dictionary){
		print "\t". $sample;
	}
print "\n";
#print data
#peakid count1 count2 ... countN
foreach my $peak_id (sort {$a<=>$b} keys %count_matrix){
	print $peak_id;
	foreach my $sample (sort keys %{ $count_matrix{$peak_id} }){
		print "\t". $count_matrix{$peak_id}{$sample};
	}
	print "\n";
}


