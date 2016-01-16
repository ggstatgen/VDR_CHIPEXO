#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#downloaded gene data from biomart like follows
#Chromosome Name Gene Start (bp) Gene End (bp)   Strand
#22      24313554        24322660        -1
#13      50202435        50208008        1
#18      2537524 2571508 -1
#20      3189514 3204516 1
#18      21178890        21242849        -1
#13      31191830        31233686        1

#I want
#chr22	24313554	24322660	+
#
my $infile;
GetOptions(
        'i=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: do_ensemblannotation_to_ucscannotation.pl -i=<FILE>\n";
     print "<FILE> input bed file\n";
     exit 1;
}

open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp $_;
	if($_ =~ /^\d{1,2}\t/){ #1,2,...,22
		my ($chr, $start, $stop, $strand) = (split /\t/)[0,1,2,3];
		$chr = 'chr' . $chr;
		if($strand eq '1'){
			$strand = '+';
		}elsif($strand eq '-1'){
			$strand = '-';
		}else{
			print "$_: strand $strand unrecognised\n";
		}
		print "$chr\t$start\t$stop\t$strand\n";
		next;
	}elsif($_ =~ /^[X|Y|M]\t/){ # X or Y
		my ($chr, $start, $stop, $strand) = (split /\t/)[0,1,2,3];
                $chr = 'chr' . $chr;
                if($strand eq '1'){
                        $strand = '+';
                }elsif($strand eq '-1'){
                        $strand = '-';
                }else{
			print "$_: strand $strand unrecognised\n";
		}
		print "$chr\t$start\t$stop\t$strand\n";
		next;
	}else{
		print $_ . "\n";
		next;
		#print "$_: unrecognised\n";
		#exit -1;
	}
}
close $instream;
