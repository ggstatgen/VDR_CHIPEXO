#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Algorithm::Combinatorics qw(combinations);

#combine all pairs of files. Send each pair to the do_macs2_idr_pair.pl script. Send each pl script to a cluster node. Use medium nodes, there are several
#the number of final jobs should be 351

my $dir;
my $SCRIPT_PATH  = '/net/isi-backup/giuseppe/scripts';
my $CODE_PATH = "/net/isi-scratch/giuseppe/tools/picard-tools-1.102";
my $qsub_opt_v = "BASH_ENV=~/.bashrc";
my $script_basename = "idr_pairs_";
GetOptions(
	'd=s'  => \$dir,
);
if(!$dir){
     print "USAGE: do_macs2_idr_pairs_cluster.pl -d=<DIR>\n";
     print "<DIR> directory containing the input .bam files\n";
     exit 1;
}


chdir $dir;
my @files = <*.bam>;
#my $combinations = combinations(\@files,2); 
my $counter = 1;

foreach my $file (@files){
	
}


while (my $pair = $combinations->next) {
	#print @$pair[0], "\t", @$pair[1], "\n";
	next if (@$pair[0] eq @$pair[1]); #I don't want to run the IDR between two instances of the same sample

	my $ID = $script_basename . $counter;
	my $outfile = $ID . rand() .  ".sh";
	my $PATH = $dir . '/' . $outfile;
	#print $PATH, "\n"; 
	
	#create bash script
	my $command = $SCRIPT_PATH . '/do_macs2_idr_pair.pl -i1=' . $dir . '/' . @$pair[0] .  ' -i2=' . $dir . '/' . @$pair[1];
	open (my $fh,  q{>}, $PATH) or die("Unable to open $PATH : $!");
	print $fh "#!/bin/bash\n";
	print $fh $command, "\n";
	close $fh;
	
	#submit bash script
	my $qsub_err = $dir . '/' . $ID . '.err';
	my $qsub_out = $dir . '/' . $ID . '.out';
	#my $submit_strg = "nice -5 qsub -v $qsub_opt_v -e $qsub_err -o  $qsub_out -q newnodes.q $PATH";
	#print $submit_strg, "\n";
	system "nice -5 qsub -v $qsub_opt_v -e $qsub_err -o  $qsub_out -q newnodes.q $PATH";
	system "rm $PATH";
	
	$counter++;
}
