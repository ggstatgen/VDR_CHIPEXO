#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Algorithm::Combinatorics qw(combinations);

#combine all pairs of files. Send each pair to the do_IDR_04_idr_pair. script. Send each pl script to a cluster node. Use medium nodes, there are several
#the number of final jobs should be 351

#the command is
#Rscript batch-consistency-analysis.r /peaks/reps/chipSampleRep1_VS_controlSampleRep0.regionPeak /peaks/reps/chipSampleRep2_VS_controlSampleRep0.regionPeak -1 /consistency/reps/chipSampleRep1_VS_chipSampleRep2 0 F signal.value

my $dir;
my $SCRIPT_PATH  = '/net/isi-backup/giuseppe/scripts';
my $R_PATH = '/net/isi-scratch/giuseppe/tools/R-3.1.0/bin/';
my $IDR_PATH = '/net/isi-scratch/giuseppe/tools/idrCode/';

my $qsub_opt_v = "BASH_ENV=~/.bashrc";
my $script_basename = "idr_pairs_";
GetOptions(
	'd=s'  => \$dir
);
if(!$dir){
     print "USAGE: do_IDR_04_idr_pairs_cluster.pl -d=<DIR>\n";
     print "<DIR> directory containing the input .regionPeak files\n";
     exit 1;
}


chdir $dir;
my @files = <*.regionPeak>;
my $combinations = combinations(\@files,2); 
my $counter = 1;
#my $chunk_counter = 400;
while (my $pair = $combinations->next) {
	#print @$pair[0], "\t", @$pair[1], "\n";
	next if (@$pair[0] eq @$pair[1]); #I don't want to run the IDR between two instances of the same sample

	my $ID = $script_basename . $counter;
	my $outfile = $ID . rand() .  ".sh";
	#my $outfile = $ID .  ".sh";
	my $PATH = $dir . '/' . $outfile;
	#print $PATH, "\n"; 
	
	#create bash script
	my $ID_PATH = $dir . '/' . $ID;
	my $command = $R_PATH . 'Rscript ' . $IDR_PATH . 'batch-consistency-analysis.r ' . $dir . '/' . @$pair[0] . ' ' .  $dir . '/' . @$pair[1] .  ' -1 ' . $ID_PATH . ' 0 F p.value';
	#aggiungi il plotting
	
	open (my $fh,  q{>}, $PATH) or die("Unable to open $PATH : $!");
	print $fh "#!/bin/bash\n";
	#print $fh "source activate\n";
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
	#if($counter % $chunk_counter == 0){
#		print "job counter: $counter\n";
#		sleep(5); #sleep half an hour
#	}
}
