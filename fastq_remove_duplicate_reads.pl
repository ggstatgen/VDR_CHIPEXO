#!/usr/bin/perl
#see http://genomics-array.blogspot.co.uk/2012/05/extracting-fastq-files-from-bamsam.html
#This is to remove redundant reads from a fastq file generated from a bam file.
#It requires that the BAM IS ORDERED so that in the fastq the redundant reads are grouped together.
use strict;
use warnings;

#my $input_file = 'test.fq';
#my $output_file = 'test_out.fq';
open (my $indata,  q{<}, $ARGV[0]) or die ("Unable to open $ARGV[0]: $!");
#open (my $outdata,  q{>}, $output_file) or die ("Unable to open $output_file: $!");

my $prev = "";

while(<$indata>){
	if(/\@(\w+):(\d+:\d+:\d+:\d+)\#\d+/){ #@DGM97JN1_120615_0211_AC0R2FACXX:3:1303:4418:86792#0
		if($prev eq $2){
        # Skips reading 3 lines
			for(0..2){
				<$indata>
    		}
    		next;
		}
        $prev = $2;
	}
	print $_;
#	print $outdata $_;
}

close $indata;
#close $outdata;
