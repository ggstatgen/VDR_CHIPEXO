#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#pre-requisite to do_impute2_cluster.pl
#script to run vcf2impute_gen and obtain a list of genotype files compatible with IMPUTE2 from a VCF. Each output will be by chromosome and computed on a different node

#options:
#vcf2impute_gen.pl
#  [-vcf]           VCF file to process; can be phased or unphased
#  [-gen]           name of new file to print; automatically gzipped
#  <-chr>           chromosome to include in output files, in (chr)[1-22,X]
#  <-start>         first position to include in output files
#  <-end>           last position to include in output files
#  <-samp_include>  file containing a list of samples to include in gen file
#  <-samp_exclude>  file containing a list of samples to exclude from gen file
#  <-samp_snptest>  flag: print sample file in SNPTEST format
#  <-special_gt>    flag: allow 'special' genotype coding in VCF
#  <-gt_like>       flag: print genotype likelihoods rather than hard calls
#  <-no_mono>       flag: omit monomorphic sites from processed files
#  <-snps_only>     flag: omit any sites that are not biallelic SNPs
#  <-indels_only>   flag: omit any sites that are not biallelic INDELs
#  <-svs_only>      flag: omit any sites that are not biallelic SVs
#
#  (args in square brackets required; args in pointy brackets optional)


my $CHR_SIZE_FILE = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes';
my $SCRIPT_PATH = '/net/isi-backup/giuseppe/scripts';
my $qsub_opt_v = "BASH_ENV=~/.bashrc";
my $script_basename = "impute2_vcf2gen_";
my $infile_vcf;

GetOptions(
        'vcf=s'         =>\$infile_vcf,
);
if(!$infile_vcf){
     print "USAGE: do_impute2_vcf_to_gen_cluster.pl -vcf=<VCF_FILE>\n";
     print "<VCF_FILE> Hapmap .vcf\n";
     exit 1;
}
my ($filename, $directory) = fileparse($infile_vcf);
$filename =~ s/(.*)\..*/$1/;



#slurp chromosome size file
my %chr_to_size;
open (my $instream,     q{<}, $CHR_SIZE_FILE) or die("Unable to open $CHR_SIZE_FILE : $!");
while(<$instream>){
	chomp;
	my ($chr, $size) = 	split(/\t/, $_);
	$chr_to_size{$chr} = $size;	
}
close $instream;

foreach my $chr (sort keys %chr_to_size){
	next if($chr =~ /random/);
	next if($chr =~ /gl/);
	next if($chr =~ /chrUn/);
	next if($chr =~ /hap/);
	
	my $size = $chr_to_size{$chr};
	my $outfile  = $directory . $script_basename . $chr . '.sh' ;
	my $qsub_err = $directory . $script_basename . $chr . '.err';
	my $qsub_out = $directory . $script_basename . $chr . '.out';
	
	#outfiles go in the infile dir
	my $outfile_gen = $directory . 'impute2_' . $filename . '_' . $chr;
	my $command = $SCRIPT_PATH . '/vcf2impute_gen -vcf ' .  $infile_vcf .  ' -gen ' .  $outfile_gen .  ' -chr ' . $chr .  ' -start 1' . ' -end ' . $size;
	open (my $fh,  q{>}, $outfile) or die("Unable to open $outfile : $!");
	print $fh "#!/bin/bash\n";
	print $fh $command, "\n";
	close $fh;

	system "nice -5 qsub -v $qsub_opt_v -e $qsub_err -o  $qsub_out -q newnodes.q $outfile";
	system "rm $outfile";
}
