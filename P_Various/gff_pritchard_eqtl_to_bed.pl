#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#convert the gff3 from here
#http://eqtl.uchicago.edu/help.html

#to bed
#very basic just keeps chr, start, end, annotation

#from gff
#chr1	Degner2012_dsQTL	Degner_dsQTL	801099	801099	3.68973163336755	.	.	DegnerDSQTL "chr1.801099 a dsQTL for 802000 to 802100" ; Note "Acts in cis "; Plot=<a href="http://eqtl.uchicago.edu/dsQTL_data/FIGURES/caQTL_7_FEB_2012/DegnerDSQTL.chr1.801099.html" >See QTL plots and view/leave comments</a>
#chr1	Degner2012_dsQTL	Degner_dsQTL	846446	846446	5.69615622511135	.	.	DegnerDSQTL "chr1.846446 a dsQTL for 846400 to 846500" ; Note "Acts in cis "; Plot=<a href="http://eqtl.uchicago.edu/dsQTL_data/FIGURES/caQTL_7_FEB_2012/DegnerDSQTL.chr1.846446.html" >See QTL plots and view/leave comments</a>
#chr1	Degner2012_dsQTL	Degner_dsQTL	901912	901912	3.94005811193805	.	.	DegnerDSQTL "chr1.901912 a dsQTL for 901300 to 901400" ; Note "Acts in cis "; Plot=<a href="http://eqtl.uchicago.edu/dsQTL_data/FIGURES/caQTL_7_FEB_2012/DegnerDSQTL.chr1.901912.html" >See QTL plots and view/leave comments</a>
#chr1	Degner2012_dsQTL	Degner_dsQTL	901458	901458	9.22112552799726	.	.	DegnerDSQTL "chr1.901458 a dsQTL for 901400 to 901500" ; Note "Acts in cis "; Plot=<a href="http://eqtl.uchicago.edu/dsQTL_data/FIGURES/caQTL_7_FEB_2012/DegnerDSQTL.chr1.901458.html" >See QTL plots and view/leave comments</a>
#chr1	Degner2012_dsQTL	Degner_dsQTL	904803	904803	4.05340037498487	.	.	DegnerDSQTL "chr1.904803 a dsQTL for 905000 to 905100" ; Note "Acts in cis "; Plot=<a href="http://eqtl.uchicago.edu/dsQTL_data/FIGURES/caQTL_7_FEB_2012/DegnerDSQTL.chr1.904803.html" >See QTL plots and view/leave comments</a>


#to bed
#chr1	1299	1300	DegnerDSQTL "chr1.801099 a dsQTL for 802000 to 802100" ; Note "Acts in cis "; Plot=<a href="http://eqtl.uchicago.edu/dsQTL_data/FIGURES/caQTL_7_FEB_2012/DegnerDSQTL.chr1.801099.html" >See QTL plots and view/leave comments</a> 5.666
#chr1	1049	1500	exon 5.6
#chr1	2999	3902	exon 4.5
#chr1	4999	5500	exon .34
#chr1	6999	9000	exon  .4345

my $infile;
GetOptions(
	'input=s'  =>\$infile,
);
if(!$infile){
     print "USAGE: gff_pritchard_eqtl_to_bed.pl -i=<INFILE>\n";
     print "<INFILE> input gff file\n";
     exit 1;
}
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
while(<$instream>){
	chomp;
	#next if ($_ =~ /##/);
	my ($chr, $start, $end, $score, $annot) = (split /\t/)[0,3,4,5,8]; 
	
	#some strange chromosome names
	next if($chr =~ /NULL/i);
	next if($chr =~ /hap/i);
	next if($chr =~ /random/i);
	if($annot =~ /(.+)\>See QTL plots and view\/leave comments\<\/a\>/){
		$annot = $1;
	}
	#need to post-process annotation. The liftover procedure does not like \" symbols etc. Replace all of them with a _
	$annot =~ s/[\s|;|"]/_/g;
	
	print $chr . "\t" . ($start-1) . "\t" . $end . "\t" . $annot . "\t" . $score . "\n";
}
close $instream;