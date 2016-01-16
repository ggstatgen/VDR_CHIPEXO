#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use List::Util qw( min max );

#13/5/2014
#carry out bayesian regression in a snp-peak association analysis setting using SNPTEST
#https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html#bayesian_tests

#things to try:
#1 different peakset
#2 different vcf set 
#3 different model
#4 multiple snps around one peak? (say +/- 200bp?)
#5 different normalisation for read count?

#Generate per-phenotype .sample file------------------------------------------------
#Obtain a sample file for SNPTEST containing the phenotypes (peak sizes)
#sample file format is detailed here:
#http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html#Sample_File_Format_

#INPUT:
#matrix of counts (obtained with do_htseq_get_matrix_from.counts.pl)
#try both normalised and not normalised
#covariate data (sex, ethnicity)

#OUTPUT:
#format follows this convention:
#ID_1 ID_2 missing cov_1 cov_2 cov_3 cov_4 pheno1 bin1
#0 0 0 D D C C P B

#with ID_1, ID_2 and missing compulsory (and corresponding 0 0 0 on second header line) followed by covariates, phenotypes, and binary
#one line per individual

#--------------------------------------
#at the moment I dynamically generate the following:
#--------------------------------------
#ID_1 ID_2 missing sex ethnicity peak_i
#0 0 0 D D P ... P
#NA06986 NA06986 0.0 male ceu 5.666
#NA06989 NA06989 0.0 male ceu 4.554
#NA06997 NA06997 0.0 female yri 2.434

#items need to be separated by SPACES, not tabs
#MISSING VALUES: code with NA

#PHENOTYPE CONTROL--------------------------------------------------------------------------
#two options:
#-quantile_normalise_phenotypes
#Quantile normalize the phenotypes. This is done AFTER samples have been excluded.
#-use_raw_phenotypes
#By default phenotypes are mean centered and scaled to have variance 1. This feature can be turned off with this option.

#INPUTS:
#see here
#https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html#input_file_formats
#genotypes can be gen, gzipped gen, bgen, vcf
#vcf has some caveats
#At the moment, I use a huge .bgen file for the genotypes (using the thresholded .gen obtained with IMPUTE2 and converted with qctools)
#when I run the code, I will use the -range option to only test the subset of SNPS under the peak in question

#summary of inputs:
#1 .bgen files with all the genotypes, split by chromosome. Provide path + basename, the program will add chromosome name
#2 .sample file with sample names, phenotypes and covariates (sex, ethnicity)

#use the script do_snptest_reads_to_sample_file.pl to generate the phenotype file
#use the .gen outputs from IMPUTE2 and the qctool to generate the .bgen genotype file

my $infile_bgen;
my $infile_rcmatrix;
my $infile_bed;
my $exclusion_list;
my $infile_pedigree = '/net/isi-scratch/giuseppe/VDR/samples.info'; #contains sex/ethnicity covariate info
GetOptions(
        'rc=s'        =>\$infile_rcmatrix,
        'b=s'         =>\$infile_bed,
		'g=s'         =>\$infile_bgen,
        'e=s'         =>\$exclusion_list
);

#$infile_bed = '/net/isi-scratch/giuseppe/VDR/POOLED_4/08_BAM_STAMPY_g1k_v37_MIN21/d_gps_smooth3_mrc20/d_bed_RBL/d_diffbind/vdr_o3_CONSENSUS.bed';
#$infile_rcmatrix = '/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_htseq_counts/d_all_30/MATRIX.txt'; #raw
#$infile_rcmatrix = '/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_htseq_counts/d_all_30/MATRIX_DESeq2_norm.txt'; #DESeq2 normalized

#Excluding Individuals (-exclude_samples,-miss_thresh)
#modify the following to decide which samples to use:
#$exclusion_list = '/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_macs2_alldup_q0.01/d_bed_RBL/d_diffbind/exclude_list.txt';
#NA10831 NA10846 NA10847 NA06997
#my $exclusion_list = '/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_htseq_counts/d_all_30/exclude_samples_29.list';
#my $exclusion_list = '/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_htseq_counts/d_all_30/exclude_samples_26.list';
#NA10846 NA11919 NA12752 NA06997
#my $exclusion_list = '/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_htseq_counts/d_all_30/exclude_samples_22.list';
#NA10846 NA11919 NA12752 NA06997 NA19214 NA19215 NA19235 NA12716

if(!$infile_rcmatrix){
     print "USAGE: do_ASSOCIATION_snptest.pl -rc=<READ_COUNT> -b=<BED_FILE> -g=<BGEN_FILE_BASENAME> -e=<EXCLUSION>\n";
     print "<READ_COUNT> matrix of read counts for the peaks in <BED_FILE> + headers. row: peak. column: sample\n";
     print "<BED_FILE> .bed file with peaks (MUST contain unique ID in NAME field (column 4)\n";
     print "(opt) <BGEN_FILE_BASENAME> path for binary .bgen genotype file (one per chr)\n";
     print "(opt) <EXCLUSION> txt file with exclusion list (one sample per line, NA format)\n";     
     print "NOTE: the genotype filename is assumed to be IMPUTE2_<chr>_subset.bgen\n";   
     exit 1;
}
if(!$infile_bed){
     print "USAGE: do_ASSOCIATION_snptest.pl -rc=<READ_COUNT> -b=<BED_FILE> -g=<BGEN_FILE_BASENAME> -e=<EXCLUSION>\n";
     print "<READ_COUNT> matrix of read counts for the peaks in <BED_FILE> + headers. row: peak. column: sample\n";
     print "<BED_FILE> .bed file with peaks (MUST contain unique ID in NAME field (column 4)\n";
     print "(opt) <BGEN_FILE_BASENAME> path for binary .bgen genotype file (one per chr)\n";
     print "(opt) <EXCLUSION> txt file with exclusion list (one sample per line, NA format)\n";     
     print "NOTE: the genotype filename is assumed to be IMPUTE2_<chr>_subset.bgen\n";   
     exit 1;
}
if(!$infile_bgen){
	$infile_bgen = '/net/isi-scratch/giuseppe/VDR/VARIANTS/IMPUTATION/d_IMPUTE_out_chksize_4mb/d_BGEN/';
	print "USING GENOTYPES FROM $infile_bgen\n\n";
}
my $TOOL_PATH       = '/net/isi-scratch/giuseppe/tools';
my $BEDTOOLS        = $TOOL_PATH . '/bedtools2-2.19.0/bin';
my $SNPTEST         = $TOOL_PATH . '/SNPTEST/snptest_v2.5_linux_x86_64_static/snptest_v2.5';
my($basename, $directory) = fileparse($infile_rcmatrix);
$basename =~ s/(.*)\..*/$1/;

my $qsub_opt_v = "BASH_ENV=~/.bashrc";
my $prefix_bgen = 'IMPUTE2_';
my $postfix_bgen = '_subset.bgen';
my $sge_job_string  = 'SNPTEST_';
my $global_out = $directory . 'snptest_' . $basename . '_global_out.txt';
my $TEMP_BASENAME = $directory . 'TEMP_SNPTEST_'; #for the temp phenotype files
my $QSUB_BASENAME = $directory . 'qsub_'; #for the sge dumps


#############
#0. slop bed intervals?
#############
#my $genome = "/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes";
#my $infile_bed_slop = $infile_bed . '.slop';
#my $SLOP = 100;
##bedtools slop [OPTIONS] -i <bed/gff/vcf> -g <genome> [-b <int> or (-l and -r)]
#system "$BEDTOOLS/bedtools slop -b $SLOP -i $infile_bed -g $genome > $infile_bed_slop";

#############
#1. get peaks and coordinates from bed file
#chr, start, stop, peakID
#this should account for any duplicates in the bed file
#############
my %peak_data;
open (my $instream,     q{<}, $infile_bed) or die("Unable to open $infile_bed : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	unless( ($_ =~ /^chr/) || ($_ =~ /^[\d+|X|Y|M]/)  ){
		print "Error in .bed file $_: either there's a header or this is not UCSC chromosome naming. Aborting..\n";
		exit;
	}
	my ($chr, $start, $stop, $ID ) = (split /\t/)[0,1,2,3];
	next if (!$chr);
	#TODO I have peaks in the following, however I don't have GTs..
	next if($chr =~ /random/);
	next if($chr =~ /gl/i);
	next if($chr =~ /Un/);
	next if($chr =~ /hap/);
	next if($chr =~ /X/);
	next if($chr =~ /Y/);
	next if($chr =~ /M/);
	
	if( ($stop - $start) <= 0 ){
		print "Error: $chr - $start - $stop is not an interval\n";
		next;
	}
	my $coord =  $start . "\t" . $stop;
	$chr = 'chr' . $chr unless ($chr =~ /^chr/); #needed for b37 data
	
	$peak_data{$chr}{$ID}{'coord'} = $coord;
}
close $instream;

#############
#2. get the matrix of read counts in the peak data structure indexed by chromosome, peak ID, sample
#############
#raw:
#peakID  NA06986 NA06989 ...
#1       1       3       ...
#normalized with DEseq
##       "NA06986"       "NA06997" 
#"1"     1.3514186257624 1.21539204543736 
#normalized with DiffBind:
##       NA06986 NA06989
#1       1.9061530108    8.6348169016
#2       1.9061530108    20.1479061037 
my @sample_names;
open ($instream,     q{<}, $infile_rcmatrix) or die("Unable to open $infile_rcmatrix : $!");
while(<$instream>){
	chomp;
	my %sample_name_to_rd;
	if($_ =~ /^\#/){
		@sample_names = split("\t", $_);
		shift @sample_names;
		#foreach (@sample_names){ s/\"//g; } # only for DEseq
		next;
	}
	my @counts = split("\t", $_);
	my $this_peak_ID = shift @counts;
	#$this_peak_ID =~ s/\"//g; # only for DEseq
	@sample_name_to_rd{@sample_names} = @counts;
	
	#need to find the chromosome where the peakID is
	foreach my $chr (keys %peak_data){
		if($peak_data{$chr}{$this_peak_ID}){
			#what about getting the maximum normalised peak rd and show it in the final file?
			$peak_data{$chr}{$this_peak_ID}{'max_rd'} = max @counts; #or maybe the median?
			foreach my $s_name (sort keys %sample_name_to_rd){
				$peak_data{$chr}{$this_peak_ID}{'sample_rd'}{$s_name} = $sample_name_to_rd{$s_name};
			}
		}
	}
}
close $instream;

################
#3 get covariate info
################
my %sample_covariates;
open ($instream,     q{<}, $infile_pedigree) or die("Unable to open $infile_pedigree : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	unless($_ =~ /^NA/){
		print "Error in covariate .info file $_: format not recognised. Must be SAMPLE,ETHNICITY,SEX - Aborting..\n";
		exit;
	}
	my ($sample_ID, $sample_ethnicity, $sample_sex ) = (split /,/)[0,1,2];
	#add to data structure
	$sample_covariates{$sample_ID}{'ethnicity'} = $sample_ethnicity;
	$sample_covariates{$sample_ID}{'sex'} = $sample_sex;	
}
close $instream;

################
#4 run genotype-phenotype association analysis
################
#you need to loop by chromosome, then by peak
#the reason why you must loop by chromosome is that the -range option only takes coordinates, not chromosome number
#so you pass it the coordinates within the chromosome
#the complication is that the genotype files must be split by chr
my %snptest_data_struct;
foreach my $peak_chr (sort keys %peak_data ){
	#build path to genotype file
	my $args_gen = $infile_bgen . $prefix_bgen . $peak_chr . $postfix_bgen;
	print "\n\nAnalysing associations for chromosome: $peak_chr\n\n";
	
	foreach my $peak_ID (sort keys %{ $peak_data{$peak_chr} }){
		#cluster-level files: three files like these per peak
		my $outfile_script  = $directory . $sge_job_string . $peak_chr . '_' . $peak_ID . '.sh' ;
		my $qsub_err        = $QSUB_BASENAME . $peak_chr . '_' . $peak_ID . '.err';
		my $qsub_out        = $QSUB_BASENAME . $peak_chr . '_' . $peak_ID . '.out';

		
		#get peak coordinates
		my ($start, $stop) = split("\t", $peak_data{$peak_chr}{$peak_ID}{'coord'});
		#build outfiles
		my $output = $directory . 'SNPTEST_' . $peak_chr . '_peak'  . $peak_ID . '.out';
		my $outlog = $directory . 'SNPTEST_' . $peak_chr . '_peak'  . $peak_ID . '.log';
		
		############################
		#4 Build per-peak phenotype file
		############################
		#--------------------------------------
		#at the moment I dynamically generate the following:
		#--------------------------------------
		#ID_1 ID_2 missing sex ethnicity peak_i
		#0 0 0 D D P ... P
		#NA06986 NA06986 0.0 male ceu 5.666
		#NA06989 NA06989 0.0 male ceu 4.554
		#NA06997 NA06997 0.0 female yri 2.434
		#generate chr-peak specific temp phenotype file
		my $temp_phenotype_file	= $TEMP_BASENAME . $peak_chr . '_peak' . $peak_ID .  'association.sample';
		open (my $temp_phenotype_stream,     q{>}, $temp_phenotype_file) or die("Unable to open $temp_phenotype_file : $!");
		my $header_1 = 'ID_1 ID_2 missing sex ethnicity pheno';
		my $header_2 = '0 0 0 D D P'; 
		print $temp_phenotype_stream $header_1, "\n";
		print $temp_phenotype_stream $header_2, "\n";
		foreach my $sample (sort keys %{ $peak_data{$peak_chr}{$peak_ID}{'sample_rd'} }){
			print $temp_phenotype_stream  $sample . ' ' . 
			                              $sample . ' 0.0 ' . 
			                              $sample_covariates{$sample}{'sex'} . ' ' . 
			                              $sample_covariates{$sample}{'ethnicity'} . ' ' . 
			                              $peak_data{$peak_chr}{$peak_ID}{'sample_rd'}{$sample} . "\n";
		}
		close $temp_phenotype_stream;
		my $args_range = $start . '-' . $stop;
		
		##########################
		#5 generate SNPTEST cluster command
		#########################
		#-quantile_normalise_phenotypes \\  
		#-exclude_samples $exclusion_list \\
		#-use_raw_phenotypes \\ 
		#-prior_qt_V_b 0.02 \\
		my $command;
		if($exclusion_list){
			$command = "$SNPTEST \\
			        -data $args_gen $temp_phenotype_file \\
			        -o $output \\
			        -bayesian 1 \\
					-range $args_range \\
			        -method expected \\
			        -cov_names ethnicity sex \\
			        -exclude_samples $exclusion_list \\
			        -pheno pheno \\
			        -quantile_normalise_phenotypes \\
			        -prior_qt_mean_b 0 \\
					-prior_qt_V_b 0.2 \\
					-prior_qt_a 3 \\
					-prior_qt_b 2  \\
			        -printids \\
			        -log  $outlog \\
			        -lower_sample_limit 10";			
		}else{
			$command = "$SNPTEST \\
		        -data $args_gen $temp_phenotype_file \\
		        -o $output \\
		        -bayesian 1 \\
				-range $args_range \\
		        -method expected \\
		        -cov_names ethnicity sex \\
		        -pheno pheno \\
		        -quantile_normalise_phenotypes \\
		        -prior_qt_mean_b 0 \\
				-prior_qt_V_b 0.2 \\
				-prior_qt_a 3 \\
				-prior_qt_b 2  \\
		        -printids \\
		        -log  $outlog \\
		        -lower_sample_limit 10";			
		}
		        

		open (my $fh,  q{>}, $outfile_script) or die("Unable to open $outfile_script : $!");
		print $fh "#!/bin/bash\n";
		print $fh $command, "\n";
		close $fh;
		
		#####################
		#6 run per-peak command on cluster
		#####################
		system "nice -5 qsub -v $qsub_opt_v -e $qsub_err -o $qsub_out -q medium_jobs.q $outfile_script";
		#system "nice -5 qsub -v $qsub_opt_v -q medium_jobs.q $outfile_script";
		unlink $outfile_script;
	}
	#ONLY when the jobs for this chromosome are complete you will delete the phenotype files
	my $COMPLETE;
	my $JOB_STATUS = `qstat`;
	print $JOB_STATUS . "\n";
	
	if ($JOB_STATUS =~ /\Q$sge_job_string/){
		$COMPLETE=0;
	}else{
		$COMPLETE=1;
	}
	while($COMPLETE == 0){
		print "CHROMOSOME $peak_chr: Waiting 10 seconds for SNPTEST jobs to complete..\n";
		sleep 10;
		$JOB_STATUS = `qstat`;
		if ($JOB_STATUS =~ /\Q$sge_job_string/){
			$COMPLETE=0;
		}else{
			$COMPLETE=1;
		}			
	}
	
	print "CHROMOSOME $peak_chr: all SNPTEST chromosome jobs completed.\n";
	system "rm $TEMP_BASENAME*.*"; #rm temp phenotypes
	system "rm $QSUB_BASENAME*.*"; #rm qsub temps
}


################
#6 gather all results strings in a unique data structure
################
#you need to open and close all files of the form SNPTEST_chr*.out 
#you need to get chr and peakID from the filename (format: SNPTEST_chr1_peak1.out)
my $sub_header; 
my $HAS_HEADER;
chdir $directory;
my @files = <SNPTEST_chr*.out>;
foreach my $file (@files) {
	my $chr = '';  my $peak_ID = '';
	if($file =~ /SNPTEST_(chr\d+)_peak(\d+)\.out/){ # only using autosomes
		$chr = $1;
		$peak_ID = $2;
	}else{
		print "ERROR: out file naming format: $file not recognised. Aborting..\n";
		exit -1;
	}
	
	open (my $outfile_stream,     q{<}, $file) or die("Unable to open $file : $!");
	while(<$outfile_stream>){
		chomp;
		next if($_ =~ /^\#/); #header 1
		if( ($_ =~ /alternate_ids rsid chromosome position alleleA alleleB index average_maximum_posterior_call info/)){  #header 2
			if(!$HAS_HEADER){
				$sub_header = $_;
				$HAS_HEADER = 1;
				next;
			}else{
				next;
			}
		}
		#discard cases where no association could be drawn (say because all samples had same genotype)
		next if($_ =~ /design_matrix_singular_value_below_limit/);
		my $data_string = $_;
		$snptest_data_struct{$chr}{$peak_ID} = $data_string;
	}
	close $outfile_stream;
}

my $final_header = 'chr peak_id max_peak_rd ' . $sub_header;
open (my $global_output_stream,     q{>}, $global_out) or die("Unable to open $global_out : $!");

print $global_output_stream $final_header, "\n";
foreach my $chr (sort keys %snptest_data_struct){
	foreach my $peak ( sort {$a<=>$b} keys %{ $snptest_data_struct{$chr} } ){
		#print $global_output_stream $chr . ' ' . $peak . ' ' . $snptest_data_struct{$chr}{$peak} . "\n";
		print $global_output_stream $chr  . ' ' . 
									$peak . ' ' .  
									$peak_data{$chr}{$peak}{'max_rd'} . ' ' .  
									$snptest_data_struct{$chr}{$peak} . "\n";
	}
}
close $global_output_stream;

#move .log and .out files to separate directories 
my $log_dir = $directory . 'd_log';
my $out_dir = $directory . 'd_out';
system "mkdir $log_dir";
system "mkdir $out_dir";
system "mv SNPTEST_*.out $out_dir";
system "mv SNPTEST_*.log $log_dir";

print "***DONE***\n";
