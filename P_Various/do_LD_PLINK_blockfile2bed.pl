#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#convert the output of plink --block command to a bed
#input sample:
# CHR          BP1          BP2           KB  NSNPS SNPS
#  22     16050678     16051882        1.205      6 rs139377059|rs6518357|rs62224610|rs79725552|rs141578542|rs2843213
#  22     16052513     16052618        0.106      2 rs149255970|rs148015672
#  22     16052838     16053435        0.598      3 rs2105538|rs147891127|rs141087528
#  22     16053509     16054800        1.292      9 rs142325435|rs141264943|rs141794080|rs11703994|rs190464175|rs139976143|rs3013006|rs55926024|rs10154476
#  22     16056965     16057288        0.324      3 rs145922255|rs192918419|rs111809146
#  22     16061100     16061294        0.195      3 rs2844848|rs140028160|rs62224636
#  22     16067555     16069707        2.153      2 rs141208996|rs1811028
#  22     16228784     16229118        0.335      2 rs142940700|rs118061773
#  22     16235732     16250859       15.128      3 rs144666531|22:16242273|rs115879525
#  22     16252594     16265053        12.46      2 rs144755802|rs75783648
#  22     16285976     16285997        0.022      2 22:16285976|rs76773885
#  22     16286465     16286468        0.004      2 rs150703810|rs73387903
#  22     16288739     16288776        0.038      3 rs78888200|rs76462367|rs67775324

#converts also 23 in X


#output: standard bed chrom start end name score
#name, the snp list
#score, the KB of the block

my $infile_plink;;
GetOptions(
		'i=s'        =>\$infile_plink,

);
if(!$infile_plink){
	print "USAGE: do_LD_plink_blockfile2bed.pl -i=<PLINK_INTERVALS>\n";
    print "<PLINK_INTERVALS> file blocks.det from plink --blocks command\n";
    exit 1;
}

open (my $instream,  q{<}, $infile_plink) or die("Unable to open $infile_plink : $!");
while(<$instream>){
	chomp;
	next if ($_ =~ /CHR/); #header
	my @fields = split(" ", $_);
	my $name =  '(' . $fields[4] . ')' . $fields[5];
	
	if($fields[0] =~ /^23/){ $fields[0] = 'X' };
	
	print $fields[0] . "\t" . $fields[1] . "\t" . $fields[2] . "\t" . $name . "\t" . $fields[3] . "\n"; 
}
close $instream;


