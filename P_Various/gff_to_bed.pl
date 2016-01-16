#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#convert gff to bed
#very basic just keeps chr, start, end, annotation

#from gff
#chr1	Canada	exon	1300	1300	.	+	.	ID=exon00001;score=1
#chr1	USA	exon	1050	1500	.	-	0	ID=exon00002;Ontology_term="GO:0046703";Ontology_term="GO:0046704"
#chr1	Canada	exon	3000	3902	.	?	2	ID=exon00003;score=4;Name=foo
#chr1	.	exon	5000	5500	.	.	.	ID=exon00004;Gap=M8 D3 M6 I1 M6
#chr1	.	exon	7000	9000	10	+	1	ID=exon00005;Dbxref="NCBI_gi:10727410"

#to bed
#chr1	1299	1300	exon 5.666
#chr1	1049	1500	exon 5.6
#chr1	2999	3902	exon 4.5
#chr1	4999	5500	exon .34
#chr1	6999	9000	exon  .4345

my $infile;
GetOptions(
	'input=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: gff_to_bed.pl -i=<INFILE>\n";
     print "<INFILE> input gff file\n";
     exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	next if ($_ =~ /##/);
	my ($chr, $annot, $start, $end, $score) = (split /\t/)[0,2,3,4,5];
	print $chr . "\t" . ($start-1) . "\t" . $end . "\t" . $annot . "\t" . $score . "\n";
}
close $instream;