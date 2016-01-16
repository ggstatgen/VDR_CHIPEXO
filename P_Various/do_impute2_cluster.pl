#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use List::Util qw(min max);
use POSIX;

#script to run impute2 and impute from the 1000g reference panel to the 14/30 samples lacking 1000g genotype
#Manchini: All the omni2.5 SNPs should be in the 1000GP dataset, so i would impute the 14 samples that have data just from HapMap, from the 1000GP haplotypes. Does that make sense? Since HapMap is phased you could do this reasonably quickly using pre-phasing

#before using this make sure you ran do_impute2_vcf2gen_cluster.pl to generate per-chromosome genotype files for the test set.

#required options
#full option list is:
#http://mathgen.stats.ox.ac.uk/impute/impute_v2.html#options

#1) -int 1 5000000 you need to split the chromosomes into "chunks". 5MB should be fine. You will need a chromosome size file for the reference you are using
#2) -g <file> this is the test data: you derive it from a .vcf file using the supplied script: vcf2impute_gen
#3) -m <file> Fine-scale recombination map for the region to be analyzed (provided in PATH_REFERENCE_DATA, genetic_map_chr[X]_combined_b37.txt)
#4) -h <file 1> <file 2> file of known haplotypes, with one row per SNP and one column per haplotype (provided in PATH_REFERENCE_DATA)
#5) -l <file 1> <file 2> Legend file(s) with information about the SNPs in the -h file(s) (provided in PATH_REFERENCE_DATA)
#6) -strand_g File showing the strand orientation of the SNP allele codings in the -g file, relative to a fixed reference point. 
#7) -Ne effective size of the population - they recommend 20000

#Example
#./impute2  
#-m ./Example/example.chr22.map 
#-h ./Example/example.chr22.1kG.haps 
#-l ./Example/example.chr22.1kG.legend 
#-g ./Example/example.chr22.study.gens 
#-strand_g ./Example/example.chr22.study.strand 
#-int 20.4e6 20.5e6 
#-Ne 20000 
#-o ./Example/example.chr22.one.phased.impute2

#il parametro int e' molto importante. Devi dividere ogni cromosoma in un intervallo di non piu di 5MB:
#The -int parameter provides an easy way to break a chromosome into smaller chunks for analysis by IMPUTE2. 
#For example, if we wanted to split a chromosome into 5-Mb regions for analysis, we could specify 
#"-int 1 5000000" for the first run of the algorithm, "-int 5000001 10000000" for the second run, and so on, all without changing the input files
#crea un ciclo esterno in cui dividi l'imputation by chr, ed un ciclo interno in cui dividi l'imputation by chunk

#once you have split a chromosome into multiple chunks and imputed them separately, the IMPUTE2 output format makes it easy to synthesize your results 
#into a single whole-chromosome file. On linux-based systems, you can simply type a command like this:
#cat chr16_chunk1.impute2 chr16_chunk2.impute2 chr16_chunk3.impute2 > chr16_chunkAll.impute2

#NOTICE: for chrX you cannot use this, you need to use other options


my $PATH_CODE = '/net/isi-scratch/giuseppe/tools/impute_v2.3.0_x86_64_static';
my $PATH_REFERENCE_DATA = '/net/isi-scratch/giuseppe/VDR/VARIANTS/IMPUTATION/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono/';
#my $PATH_REFERENCE_DATA = '/net/isi-scratch/giuseppe/VDR/VARIANTS/IMPUTATION/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing/'; #quale usare? vedi
my $CHR_SIZE_FILE = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes';

#hardcoded impute reference data basenames
#map: format is genetic_map_chr10_combined_b37.txt
#need pre and post fix
my $MAP_FILE_PRE  =  'genetic_map_';
my $MAP_FILE_POST =  '_combined_b37.txt'; 
#haplo: format is ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz
#need pre and post fix
my $REFDATA_FILE_PRE = 'ALL.';
my $HAPLO_FILE_POST = '.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz';
#legend: format is ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz
my $LEG_FILE_POST = '.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz';

my $qsub_opt_v = "BASH_ENV=~/.bashrc";

my $infile_gen;
my $chunksize;
my $start;
my $stop;
GetOptions(
        'gen=s'         =>\$infile_gen,
        'chunksize=i'   =>\$chunksize
);
if(!$infile_gen){
     print "USAGE: do_imputer2_cluster.pl -gen=<GEN_FILE_BASENAME> -chunksize=<CHUNKSIZE>\n";
     print "<GEN_FILE> path + basename for gzipped Hapmap test data file to impute. Obtained from .vcf via vcf2impute_gen. eg 'impute2_genotypes_CEUYRI_30_r27_nr.hg19_fwd'\n";
     print "<CHUNKSIZE> impute2 chunk in bytes. 5Mb or shorter recommended (Eg 5000000, 3000000, etc)\n";
     print "NOTE: script expects same basename for each chromosome file. It will add '_chrX.gz'\n";
     exit 1;
}
if(!$chunksize){
     print "USAGE: do_imputer2_cluster.pl -gen=<GEN_FILE_BASENAME> -chunksize=<CHUNKSIZE>\n";
     print "<GEN_FILE> path + basename for gzipped Hapmap test data file to impute. Obtained from .vcf via vcf2impute_gen. eg 'impute2_genotypes_CEUYRI_30_r27_nr.hg19_fwd'\n";
     print "<CHUNKSIZE> impute2 chunk in bytes. 5Mb or shorter recommended (Eg 5000000, 3000000, etc)\n";
     print "NOTE: script expects same basename for each chromosome file. It will add '_chrX.gz'\n";
     exit 1;
}

my ($filename, $directory) = fileparse($infile_gen);
#slurp chromosome_size file
my %chr_to_size;
open (my $instream,     q{<}, $CHR_SIZE_FILE) or die("Unable to open $CHR_SIZE_FILE : $!");
while(<$instream>){
	chomp;
	my ($chr, $size) = 	split(/\t/, $_);
	$chr_to_size{$chr} = $size;	
}
close $instream;

#---------------------
#LOOP 1 -  chromosome level
#---------------------
#chrX requires a special command and a special sample file to tell females from males
my $SAMPLEG_FILE = $directory . 'impute2_genotypes_CEUYRI_30_r27_nr.hg19_fwd.sampledata';

foreach my $chr (sort keys %chr_to_size){
	next if($chr =~ /random/);
	next if($chr =~ /gl/);
	next if($chr =~ /chrUn/);
	next if($chr =~ /hap/);
	next if($chr =~ /chrX/);
	next if($chr =~ /chrY/);
	next if($chr =~ /chrM/);

	my $chr_size = $chr_to_size{$chr};
	my $this_chr_chunk_number = ceil($chr_size / $chunksize);
	my $chunk_counter = 1;
	
	#build impute2 input filenames for this chr
	#need to build based on a hardcoded basename
	my $M_FILE = $PATH_REFERENCE_DATA . $MAP_FILE_PRE     . $chr . $MAP_FILE_POST;
	my $H_FILE = $PATH_REFERENCE_DATA . $REFDATA_FILE_PRE . $chr . $HAPLO_FILE_POST;
	my $L_FILE = $PATH_REFERENCE_DATA . $REFDATA_FILE_PRE . $chr . $LEG_FILE_POST;
	my $G_FILE = $infile_gen . '_' . $chr . '.gz';
	my $STRAND_FILE; #?? dove lo prendo questo	
	
	#initialise chunk---
	my $start_chunk_coord = 1;
	my $end_chunk_coord = $chunksize;
	
	#wildcards (for log file manipulations)
	#basenames 
	#IMPUTE2_chr1_chunk75.out
	#IMPUTE2_chr1_chunk*.out_info
	#IMPUTE2_chr1_chunk*.out_info_by_sample
	#IMPUTE2_chr1_chunk*.out_summary
	my $chunkout_wildcard     = 'IMPUTE2_' . $chr . '_chunk*.gen';
	my $chunkinfo1_wildcard   = 'IMPUTE2_' . $chr . '_chunk*.gen_info';;
	my $chunkinfo2_wildcard   = 'IMPUTE2_' . $chr . '_chunk*.gen_info_by_sample';
	my $chunksummary_wildcard = 'IMPUTE2_' . $chr . '_chunk*.gen_summary';
	my $LOG_PATH = $directory . 'd_IMPUTE2_' . $chr . '/';
	
	#-------------------
	#LOOP 2 -  chunk level
	#-------------------
	my $concatenated_out_string = 'cat'; #let's concatenate by chromosome and remove the single chunk outputs
	while($chunk_counter <= $this_chr_chunk_number){
		my $command;
		my $script_basename = 'IMPUTE2_' . $chr . '_chunk' . $chunk_counter;
		#outputs
		my $outfile_script  = $directory . $script_basename . '.sh' ;
		my $qsub_err        = $directory . $script_basename . '.err';
		my $qsub_out        = $directory . $script_basename . '.out';
		my $OUT_FILE        = $directory . $script_basename . '.gen';
				
		#build the chunk interval
		my $INT = $start_chunk_coord . ' ' . $end_chunk_coord;
		
		#create IMPUTE2 command
		#if chrX, create specific command
		if($chr =~ /X/i){
			$command = "$PATH_CODE/impute2          \\
							-chrX                   \\
							-m  $M_FILE	            \\
							-h  $H_FILE	            \\
							-l  $L_FILE	            \\
							-g  $G_FILE	            \\
							-sample_g $SAMPLEG_FILE \\
 							-int $INT               \\
							-Ne 20000		        \\
							-o $OUT_FILE";			
		}else{
			$command = "$PATH_CODE/impute2          \\
							-m  $M_FILE	            \\
							-h  $H_FILE	            \\
							-l  $L_FILE	            \\
							-g  $G_FILE	            \\
							-sample_g $SAMPLEG_FILE \\
							-int $INT               \\
							-Ne 20000		        \\
							-o $OUT_FILE";			
		}

		#send this to the cluster						
		open (my $fh,  q{>}, $outfile_script) or die("Unable to open $outfile_script : $!");
		print $fh "#!/bin/bash\n";
		print $fh $command, "\n";
		close $fh;
		
		#submit bash script
		#my $submit_strg = "nice -5 qsub -v $qsub_opt_v -e $qsub_err  -q newnodes.q $outfile_script";
		#print $submit_strg, "\n";
		system "nice -5 qsub -v $qsub_opt_v -e $qsub_err -o $qsub_out -q medium_jobs.q $outfile_script";
		system "rm $outfile_script";
		
		#update chunk window
		$start_chunk_coord += $chunksize;
		$end_chunk_coord += $chunksize;
		$end_chunk_coord = min($end_chunk_coord, $chr_size); # the last chunk might be smaller than the $chunksize window
		
		$chunk_counter++;
	}
	
	#check jobs are finished for this chr, then concatenate
	my $COMPLETE;
	my $JOB_STATUS = `qstat`;
	print $JOB_STATUS . "\n";
	
	if ($JOB_STATUS =~ /IMPUTE2/){
		$COMPLETE=0;
	}else{
		$COMPLETE=1;
	}
	while($COMPLETE == 0){
		print "CHROMOSOME $chr: Waiting 60 seconds for IMPUTE2 chunk jobs to complete..\n";
		sleep 60;
		$JOB_STATUS = `qstat`;
		if ($JOB_STATUS =~ /IMPUTE2/){
			$COMPLETE=0;
		}else{
			$COMPLETE=1;
		}			
	}
	print "CHROMOSOME $chr: all IMPUTE2 chunk jobs completed.\n";
	#cat chunks before moving to the next chr
	#some chunks will be empty - build concatenate-string hash excluding the empty ones
	$chunk_counter = 1;
	while($chunk_counter <= $this_chr_chunk_number){
		my $script_basename = 'IMPUTE2_' . $chr . '_chunk' . $chunk_counter;
		my $out_filename    = $directory . $script_basename . '.gen';
		if(-e $out_filename){
			$concatenated_out_string .= ' ' . $out_filename;
		}
		$chunk_counter++;
	}
	my $this_chr_global_out = $directory . 'IMPUTE2_' . $chr . '.gen';
	system "$concatenated_out_string > $this_chr_global_out";

	#-------------------------------------
	#remove chunk outputs, move chunk logs
	#-------------------------------------
	system "rm $chunkout_wildcard";
	#move summary, info and info by sample  files in IMPUTE_chr directory
	system "mkdir $LOG_PATH";
	system "mv $chunkinfo1_wildcard $LOG_PATH";
	system "mv $chunkinfo2_wildcard $LOG_PATH";
	system "mv $chunksummary_wildcard $LOG_PATH";
}