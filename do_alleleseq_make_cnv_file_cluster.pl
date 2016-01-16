#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Algorithm::Combinatorics qw(combinations);

#This calls do_alleleseq_make_cnv_file.pl and runs it on all pairs vcf / cnv in the directory

my $dir;
my $SCRIPT_PATH  = '/net/isi-backup/giuseppe/scripts';
#my $qsub_opt_v = "BASH_ENV=~/.bashrc";
my $script_basename = "alleleseq_make_cnv_";
GetOptions(
	'd=s'  => \$dir
);
if(!$dir){
     print "USAGE: do_alleleseq_make_cnv_file_cluster.pl -d=<DIR>\n";
     print "<DIR> directory containing the input .vcf.gz file and .cnv file, for each sample\n";
     exit 1;
}


chdir $dir;
my @files = <*.*>;
my $combinations = combinations(\@files,2); 
while (my $pair = $combinations->next) {
	my $ID1; my $ID2;
	
	if(@$pair[0] =~ /.*(NA\d{5})\.vcf\.gz/){
		$ID1 = $1;
	}elsif(@$pair[0] =~ /.*(NA\d{5}).*\.tsv/){
		$ID1 = $1;
	}else{
		print "Error: file format not recognised, aborting.\n";
		exit -1;
	}
	if(@$pair[1] =~ /.*(NA\d{5}).*\.tsv/){
                $ID2 = $1;
        }elsif(@$pair[1] =~ /.*(NA\d{5})\.vcf\.gz/){
		$ID2 = $1;
	}else{
                print "Error: file format not recognised, aborting.\n";
                exit -1;
        }
	next unless ($ID1 eq $ID2);
	#iprint @$pair[0], "\t", @$pair[1], "\n";
	my $arg_vcf; my $arg_tsv;
	if(@$pair[1] =~ /.*\.tsv/){
		$arg_vcf = @$pair[0]; $arg_tsv = @$pair[1];
	}else{
		$arg_vcf = @$pair[1]; $arg_tsv = @$pair[0];
	}

	my $outfile = $script_basename . $ID1 .  ".sh";
	my $PATH = $dir . '/' . $outfile;
	
	#create bash script
	my $command = $SCRIPT_PATH . '/do_alleleseq_make_cnv_file.pl -cnv=' . $dir . '/' . $arg_tsv .  ' -vcf=' . $dir . '/' . $arg_vcf;
	open (my $fh,  q{>}, $PATH) or die("Unable to open $PATH : $!");
	print $fh "#!/bin/bash\n";
        print $fh "source activate\n";
	print $fh $command, "\n";
	close $fh;
	
	#submit bash script
	my $qsub_err = $dir . '/' . $ID1 . '.err';
	my $qsub_out = $dir . '/' . $ID1 . '.cnv';
	#my $submit_strg = "nice -5 qsub -v $qsub_opt_v -e $qsub_err -o  $qsub_out -q newnodes.q $PATH";
	#print $submit_strg, "\n";
	system "nice -5 qsub -e $qsub_err -o  $qsub_out -q medium_jobs.q $PATH";
	system "rm $PATH";
}
