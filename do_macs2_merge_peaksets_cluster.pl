#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Algorithm::Combinatorics qw(combinations);

#31/7/2014
#runs do_macs2_merge_peaksets.pl across all samples
#it should spawn 30 jobs
#

my $dir;
my $SCRIPT_PATH  = '/net/isi-backup/giuseppe/scripts';
#my $qsub_opt_v = "BASH_ENV=~/.bashrc";
my $script_basename = "pkst_mrg_";
GetOptions(
	'd=s'  => \$dir
);
if(!$dir){
     print "USAGE: do_macs2_merge_peaksets_cluster.pl -d=<DIR>\n";
     print "<DIR> directory containing the input .bed files, TWO PER SAMPLE (different d_size each)\n";
     exit 1;
}


chdir $dir;
my @files = <*.bed>;
my $combinations = combinations(\@files,2); 
while (my $pair = $combinations->next) {
	my $ID1; my $ID2;
	my $dist1; my $dist2;


	if(@$pair[0] =~ /.*(GM\d{5}).*(d\d+).*/){
		$ID1 = $1;
		$dist1 = $2;
	}else{
		print "Error: file format not recognised, aborting.\n";
		exit -1;
	}
	if(@$pair[1] =~ /.*(GM\d{5}).*(d\d+).*/){
                $ID2 = $1;
		$dist2 = $2;
        }else{
                print "Error: file format not recognised, aborting.\n";
                exit -1;
        }
	next unless ($ID1 eq $ID2);
	next if ($dist1 eq $dist2);
	#print @$pair[0], "\t", @$pair[1], "\n";

	my $outfile = $script_basename . $ID1 .  ".sh";
	my $PATH = $dir . '/' . $outfile;
	
	#create bash script
	my $command = $SCRIPT_PATH . '/do_macs2_merge_peaksets.pl -i1=' . $dir . '/' . @$pair[1] .  ' -i2=' . $dir . '/' . @$pair[0];
	open (my $fh,  q{>}, $PATH) or die("Unable to open $PATH : $!");
	print $fh "#!/bin/bash\n";
        print $fh "source activate\n";
	print $fh $command, "\n";
	close $fh;
	
	#submit bash script
	my $qsub_err = $dir . '/' . $ID1 . '.err';
	my $qsub_out = $dir . '/' . $ID1 . '_peaks_consensus.bed';
	#my $submit_strg = "nice -5 qsub -v $qsub_opt_v -e $qsub_err -o  $qsub_out -q newnodes.q $PATH";
	#print $submit_strg, "\n";
	system "nice -5 qsub -e $qsub_err -o  $qsub_out -q newnodes.q $PATH";
	system "rm $PATH";
}
