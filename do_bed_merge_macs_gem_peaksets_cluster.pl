#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Algorithm::Combinatorics qw(combinations);

#31/7/2014
#runs do_bed_merge_macs_gem_peaksets.pl across all samples
#it should spawn 30 jobs

#USAGE: do_bed_merge_macs_gem_peaksets.pl -i1=<INFILE1> -i2=<INFILE2> -d=<DISTANCE>
#<INFILE1> bed file 1
#<INFILE2> bed file 2
#(optional)<DISTANCE> whether to merge features distant <DISTANCE>bp apart. Default: do not merge non-overlapping features

#expected file name for the pair is 
#1 VDR_GM19249_peaks_macs.bed 
#2 VDR_GM19249_peaks_gem.bed

my $dir;
my $distance;
my $SCRIPT_PATH  = '/net/isi-backup/giuseppe/scripts';
my $script_basename = "pkst_mrg_";
GetOptions(
	'dir=s'  => \$dir,
	'd=i'    => \$distance
);
if(!$dir){
     print "USAGE: do_bed_merge_macs_gem_peaksets.pl -dir=<DIR> -d=<DISTANCE>\n";
     print "<DIR> directory containing the input .bed files, TWO PER SAMPLE (..GMxxxxx_peaks_macs.. and ..GMxxxxx_peaks_gem..)\n";	
     print "(optional)<DISTANCE> whether to merge features distant <DISTANCE>bp apart. Default: do not merge non-overlapping features\n";
     exit 1;
}

chdir $dir;
my @files = <*.bed>;
my $combinations = combinations(\@files,2); 
while (my $pair = $combinations->next) {
	my $ID0; my $ID1;
	my $pcaller0; my $pcaller1;

	if( (@$pair[0] =~ /.*(GM\d{5}).*(macs).*/) || (@$pair[0] =~ /.*(GM\d{5}).*(gem).*/)  ){
		$ID0 = $1;
		$pcaller0 = $2;
	}else{
		print "Error: file format not recognised, aborting.\n";
		exit -1;
	}
        if( (@$pair[1] =~ /.*(GM\d{5}).*(macs).*/) || (@$pair[1] =~ /.*(GM\d{5}).*(gem).*/)  ){
                $ID1 = $1;
		$pcaller1 = $2;
        }else{
                print "Error: file format not recognised, aborting.\n";
                exit -1;
        }
	next unless ($ID0 eq $ID1);
	next if ($pcaller0 eq $pcaller1);
	#print @$pair[0], "\t", @$pair[1], "\n";

	my $outfile = $script_basename . $ID0 .  ".sh";
	my $PATH = $dir . '/' . $outfile;
	
	#create bash script
	my $command;
	if($distance){
	        $command = $SCRIPT_PATH . '/do_bed_merge_macs_gem_peaksets.pl -i1=' . $dir . '/' . @$pair[0] .  ' -i2=' . $dir . '/' . @$pair[1] . ' -d=' . $distance;
	}else{
        	$command = $SCRIPT_PATH . '/do_bed_merge_macs_gem_peaksets.pl -i1=' . $dir . '/' . @$pair[0] .  ' -i2=' . $dir . '/' . @$pair[1];
	}

	open (my $fh,  q{>}, $PATH) or die("Unable to open $PATH : $!");
	print $fh "#!/bin/bash\n";
        print $fh "source activate\n";
	print $fh $command, "\n";
	close $fh;
	
	#submit bash script
	my $qsub_err = $dir . '/' . $ID0 . '.err';
	my $qsub_out = $dir . '/' . $ID1 . '_peaks_consensus.bed';
	#my $submit_strg = "nice -5 qsub -v $qsub_opt_v -e $qsub_err -o  $qsub_out -q newnodes.q $PATH";
	#print $submit_strg, "\n";
	system "nice -5 qsub -e $qsub_err -o  $qsub_out -q newnodes.q $PATH";
	system "rm $PATH";
}
