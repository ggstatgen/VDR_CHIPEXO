#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;
use List::Util qw(sum);
use Statistics::R;

########
#VET - variant enrichment tester - cluster version
########


#just like the other similarly named ones, but the enrichment is done by the cluster nodes, one per phenotype, calling the routine
#do_funseq_VET_02_freqmatch_bootstrap_mt_inner_noBT.pl

#SNPsnap returns 1000 sets of random snps which are frequency matched with the input
#However these sets are not under ChIPexo VDR binding coverage.
#So while they're useful because I can correct for DAF frequency matching, gene features etc, I think the produced enrichments will still not be ideal
#because I need to remove VDR binding bias.
#My foreground is biased by where my read coverage is. I had found already that VDR peaks were enriched in GWAS regions. So when looking at VDR variants, I need to factor out this enrichment by only considering variants which were tested, so variants under >5 read depth coverage.

#In this script I do all by myself. Given an LD-clumped foreground, and an [LD-clumped] background, I will generate 1000 subsets from the background which are DAF matched to the foreground. Then I will ran the enrichment tests as done in do_funseq_GAT_SNPtest_randomizer.pl
#To do this I will use a binning approach, as follows:
#http://www.cureffi.org/2013/01/29/sampling-a-matching-distribution-for-bootstrapping/
#To handle this problem, Nicolae drew random sets of SNPs matching on MAF.  Specifically, he binned all of the SNPs by increments of 0.05, so there are the 0.00 - 0.05 MAF SNPs, the 0.05 - 0.10 MAF SNPs, etc.  Then he drew the same number from each bin in the random sets as were present in the true GWAS set.

#INPUTS:
#1 bed.gz of FOREGROUND snps, LD-clumped [D' definition, CEU] using do_funseq_GAT_LDclump_FG_and_BG.pl
#2 bed.gz of BACKGROUND snps, [LD clumped? TEST this] with 1kG allele frequencies to derive DAF for all foreground and background SNPs
#3 ANNOTATION: GRASP LD blocks [D' definition, CEU]

#METHOD:
#-get DAF for all fg snps. Create a BINNING hash with twenty keys. bin1 [0.00-0.05], bin2 ]0.05-1.00], ..., ]0.95,1.00].
#gather the foreground snps in each bin based on DAF (MAF?) and count them
#compute foreground intersection with phenotypes
#then draw the same number of snps from the Background set so that THE SAME binning happens. Save the 1000 sets in a file or in hash
#for each phenotype, run enrichment test. Save stats. Repeat.
#correct pvalues

#RANDOM SELECTION OF n ELEMENTS FROM BIN:
#Using Randall Swartz solution:
#http://www.perlmonks.org/?node=randomly%20choosing%20elements%20from%20an%20array
#Shuffling is still more work than you need to do if you can be destructive to the list. To pick $n items from list @list, use:
#my @result;
#push @result, splice @list, rand @list, 1 while @result < $n;

my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";
my $PERL = "/net/isi-cgat/ifs/apps/apps/perl-5.16.1/bin/perl";
my $ENRICHMENT_ROUTINE = "/net/isi-backup/giuseppe/scripts/P_Various/do_funseq_VET_inner_noBT.pl";
my $qsub_opt_v = "BASH_ENV=~/.bashrc";
my $sge_job_string  = 'VET_JOBS_';


my $LD_PLINK_GWAScat = "/net/isi-scratch/giuseppe/indexes/GWAS_Gwascat/UCSC_gwas_catalog/LD_PLINK_CEU_GWAScat.maxp_5E-8.minsnp5.bed.gz";
#my $LD_PLINK_GWAScat = "/net/isi-scratch/giuseppe/indexes/GWAS_Gwascat/UCSC_gwas_catalog/LD_PLINK_CEU_GWAScat.maxp_5E-8.bed.gz";
my $LD_PLINK_GRASP   = "/net/isi-scratch/giuseppe/indexes/GWAS_GRASP/LD_PLINK_CEU_GRASP2_plus_Beecham2013.bed.gz";
my $LD_r2_MIG        = "/net/isi-scratch/giuseppe/indexes/GWAS_GRASP/LD_MIG_r0.8_EUR_noindels_maf0.01_hwe0.001_GRASP.bed.gz";

my $BIN_NUMBER = 20;
my $RANDOM_SETS;
my $INPUT_GWAS_BLOCKS;
#############
#globals
#############
my %phenotype_map;
my %variants_1kg;
my %DAF_binning_FG;
my %DAF_binning_BG;
my $ld_block_type;

#hash ethnicity ids
my $AFEUR  = 'AFEUR';
my $AFGLOB = 'AFGLOB';
my $AFAFR  = 'AFAFR'; 

#get the MAF/DAF values for foreground and background
#background: testing all the ~116,000 variants atm
#blocks: testing the usual D' CEU blocks on GRASP GWAS snps.
#also testing r2 MIG ld blocks at r2 > 0.8

my $INPUT_FOREGROUND;
my $INPUT_BACKGROUND;
my $ALLELE_FREQ;
my $use_gwascat;
GetOptions(
        'fg=s'      =>\$INPUT_FOREGROUND,
        'bg=s'      =>\$INPUT_BACKGROUND,
        'r=s'       =>\$RANDOM_SETS,
        'af=s'      =>\$ALLELE_FREQ,
        'ld=s'      =>\$ld_block_type,
        'gwascat'   =>\$use_gwascat
);

#$INPUT_FOREGROUND = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/ET-VDRBV-REP/LDCLUMP_REP_FG_VDRBVs_Output_noDBRECUR.bed.gz";
#$INPUT_BACKGROUND = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/d_GAT_BACKGROUNDS/GAT_BACKGROUND_SIMASYM_FINAL2.vcf.gz";
#$ALLELE_FREQ = 'maf';
#$ld_block_type = 'D';

my $USAGE = "USAGE: do_funseq_VET_02_mt.pl -fg=<FG> -bg=<BG> -r=<RANDOMSETS> -af=<daf|maf>  -ld=<D|R> -gwascat\n" . 
    		"<FG> b37 bed.gz containing the LD-CLUMPED foreground SNPs [generated by do_funseq_ET_01_LDclump_FG_and_BG.pl]\n" .
    		"<BG> b37 vcf.gz containing the [LD-CLUMPED] background SNPs with MAF information [generated by do_funseq_ET_01_LDclump_FG_and_BG.pl]\n" .
    		"<RANDOMSETS> number of random sets to gather from the background (eg 1000)\n" . 
			"<af> allele frequency to consider: one of maf (alternate allele freq) or daf (derived allele freq)\n" .
			"<ld> whether to use plink D'[D] or Mig Rsquare >0.8 [R] definition of GWAS LD blocks\n" .
			"<gwascat> if <ld> = D, then setting this will use tag snps mapped to LD blocks taken from the GWAS catalog, rather than the GRASP\n";  


if(!$INPUT_FOREGROUND){
	print $USAGE;
    exit 1;
}
if(!$INPUT_BACKGROUND){
	print $USAGE;
    exit 1;
}
if(!$RANDOM_SETS){
	print $USAGE;
    exit 1;
}
if(!$ALLELE_FREQ){
	print $USAGE;
    exit 1;
}
if(!$ld_block_type){
	print $USAGE;
    exit 1;
}
unless($ALLELE_FREQ eq 'maf' || $ALLELE_FREQ eq 'daf'){
	print $USAGE;
    exit 1;
}
unless($ld_block_type eq 'D' || $ld_block_type eq 'R'){
	print $USAGE;
    exit 1;	
}
if($use_gwascat){
	unless($ld_block_type eq 'D'){
		print $USAGE;
    	exit 1;			
	}
}

my $ld_label;
if($ld_block_type eq 'D'){
	if($use_gwascat){
		$INPUT_GWAS_BLOCKS = $LD_PLINK_GWAScat;
		$ld_label = '_LD_Dprime_GWAScat';		
	}else{
		$INPUT_GWAS_BLOCKS = $LD_PLINK_GRASP;
		$ld_label = '_LD_Dprime_GRASP';
	}
}else{
	$INPUT_GWAS_BLOCKS = $LD_r2_MIG;
	$ld_label = '_LD_r2_0.8';	
}

#difference in blocks
#---------------------D'
#1 1
#2 802496
#3 808631
#4 (9)rs10157494|rs139867617|rs142367470|rs7526310|rs72631880|rs61768215|rs4951933|rs113532475|rs11240779
#5 6.136
#6 1
#7 807512
#8 rs10751454
#9 -
#10 -
#11 -
#12 Cardiovascular disease prevalence
#13 Cardiovascular disease (CVD);Myocardial infarction (MI);Mortality;Cancer(p=2.20E-13)
#14 1

#---------------------r2
#1 1
#2 807511
#3 807512
#4 (1)rs10751454
#5 -
#6 1
#7 807512
#8 rs10751454
#9 -
#10 -
#11 -
#12 Cardiovascular disease prevalence
#13 Cardiovascular disease (CVD);Myocardial infarction (MI);Mortality;Cancer(p=2.20E-13)
#14 0


my($basename, $directory) = fileparse($INPUT_FOREGROUND);
$basename =~ s/(.*)\..*/$1/;
$basename =~ s/(.*)\..*/$1/;

my $out_results = $directory . 'RESULTS_VET_' . $basename . $ld_label . '_match_' . $ALLELE_FREQ . '_r' . $RANDOM_SETS . '.tsv';
my $out_temp_buffer         =  $directory     . $basename . '_TEMP_diseaseblocks.tsv';
my $out_temp_buffer_ph_list =  $directory     . $basename . '_TEMP_diseaseblocks_phlist.txt';
my $temp_buffer_ph_map      =  $directory     . $basename . '_TEMP_diseaseblocks_phmap.txt';

#collect phenotypes which intersect foreground snps as these will be the ones to test.

#####################
#1 -build  hash to map chr-pos to DAF-CEU and DAF-YRI (the latter not needed now)
#####################
#before saving the frequencies, you need to to KNOW if the ancestral allele is the ref or the alt
#if the ancestral allele is the ref, save the frequency as is
#if the ancestral allele is the alt, the frequency you have is for the ancestral. The derived will be (1 - freq)
#however if you choose 'maf' below, you'll simply get what's reported in the AF field, without caring about ancestral allele

if($ALLELE_FREQ eq 'daf'){
	print STDERR "Collecting 1kg DAF..\n";
	get_DAF(\%variants_1kg);	
}elsif($ALLELE_FREQ eq 'maf'){
	print STDERR "Collecting 1kg MAF..\n";
	get_MAF(\%variants_1kg);	
}else{
	print STDERR "Should not be here. Aborting..\n";
	exit -1;
}

###########
#2 intersect FG with disease LD blocks to see which traits are represented.
###########
#intersect disease blocks with foreground SNPs
system "$BEDTOOLS intersect -wo -a $INPUT_GWAS_BLOCKS -b $INPUT_FOREGROUND > $out_temp_buffer";
#collect intersecting phenotypes
system "cut -f 12 $out_temp_buffer | sort | uniq > $out_temp_buffer_ph_list";
#then get phenotype list from file
open (my $instream,  q{<}, $out_temp_buffer_ph_list) or die("Unable to open $out_temp_buffer_ph_list : $!");
while(<$instream>){
	chomp;
	next if($_ eq '');
	$phenotype_map{lc($_)} = 1;
}
close $instream;
#now we have a list of phenotypes seen in correspondence of foreground variants. Create a temp bed for each term and intersect with foreground and with each dynamically generated background bed. Store number of intersections. Calculate z-scores and fdr correct

##############
#3 - bin foreground SNPs by DAF/MAF
##############
#open the foreground snps and query them to %variants_1kg
print STDERR "Binning FG variants..\n";

my ($total_fg_vars, $skipped_fg_vars, $valid_fg_vars)  = bin_variants($INPUT_FOREGROUND,\%DAF_binning_FG, 'fg');

print STDERR "$total_fg_vars total fg vars\n";
print STDERR "$skipped_fg_vars fg vars skipped (no key found/no DAF info found)\n";
print STDERR "Valid fg variants: $valid_fg_vars\n";
print STDERR "FOREGROUND DAF HISTOGRAM:\n";

print_histogram(\%DAF_binning_FG, 'fg');

##############
#4 - bin background SNPs by DAF/MAF
##############
print STDERR "\n\nBinning BG variants..\n";

my ($total_bg_vars, $skipped_bg_vars, $valid_bg_vars)  = bin_variants($INPUT_BACKGROUND,\%DAF_binning_BG, 'bg');

print STDERR "$total_bg_vars total bg vars\n";
print STDERR "$skipped_bg_vars bg vars skipped (no key found/no DAF info found)\n";
print STDERR "Valid bg variants: $valid_bg_vars\n";
print STDERR "BACKGROUND DAF HISTOGRAM:\n";

print_histogram(\%DAF_binning_BG, 'bg');

#cleanup
%variants_1kg = ();

############
#5 - collect frequency matched variants and save them into a matrix file
############
print STDERR "Collecting $RANDOM_SETS sets of $valid_fg_vars frequency-matched background variants..\n";

#create array of array containing 1000 random frequency matched sets
#for each foreground bin, count the number of instances in the bin, then get a random set of n instances from the background bin
#if the foreground bin is empty, keep the background bin empty.

my $out_temp_bg_matrix =  $directory . $basename . '_TEMP_bgmatrix.txt';
open (my $tempfile,  q{>}, $out_temp_bg_matrix) or die("Unable to open $out_temp_bg_matrix : $!");

while ($RANDOM_SETS > 0) {
	my @this_random_set;
	
	for (0..$BIN_NUMBER){
		$_ = $_/$BIN_NUMBER;	 
	     if($DAF_binning_FG{$_}){
	     	my @all_bg_vars_this_bin;
	    	my @match_bg_vars_this_bin;  	
	    	my $n = $DAF_binning_FG{$_};
			#get n items at random from $DAF_binning_BG{$_}
			#you should chech that the background contains MORE variants than the fg
			#though you know it does for this particular bg/fg combination
			#my $n_bg = scalar(keys %{$DAF_binning_BG{$_}});
			#if($n > $n_bg).. 
			@all_bg_vars_this_bin = keys %{$DAF_binning_BG{$_}};
			#is there a way to get n random variants?
			push @match_bg_vars_this_bin, splice @all_bg_vars_this_bin, rand @all_bg_vars_this_bin, 1 while @match_bg_vars_this_bin < $n;			
			push(@this_random_set, @match_bg_vars_this_bin); #DOUBLE CHECK
	    }else{
	    	#there are no fg variants in this bin, there should be none in this bg then
	    	next;
	    }	
	}
	#print the random set
	print $tempfile join("\t",@this_random_set),"\n";
	#push @random_sample_matrix, [ @this_random_set ];
	$RANDOM_SETS -= 1;
}
close $tempfile;

#####################
#6 - create file with id mapping to observed phenotypes
#it seems some trait names are problematic. I map them to IDs and pass the ID and the map to the routine below
######################
open (my $tempstream,  q{>}, $temp_buffer_ph_map) or die("Unable to open $temp_buffer_ph_map : $!");
my $job_counter = 1;
foreach my $item (sort keys %phenotype_map){
	print $tempstream $job_counter . "\t" . $item . "\n";
	$job_counter++;
}
close $tempstream;


##########=#########
#6- COMPUTE ENRICHMENTS
####################
$job_counter = 1;
#my $temp_counter = 3;
my $QSUB_BASENAME = $directory . $sge_job_string. 'qsub_'; #for the sge dumps
foreach my $item (sort keys %phenotype_map){
	#cluster-level files: three files like these per peak
	
	my $sge_script      = $QSUB_BASENAME               . 'phenotype_' . $job_counter .  '.sh' ;
	my $qsub_err        = $QSUB_BASENAME               . 'phenotype_' . $job_counter  . '.err';
	my $qsub_out        = $QSUB_BASENAME               . 'phenotype_' . $job_counter  . '.out';

#				        -ph=\"$item\"             \\
	my $command = "source activate\n" . 
					    "$PERL $ENRICHMENT_ROUTINE     \\
				        -fg=$INPUT_FOREGROUND          \\
				        -bg=$out_temp_bg_matrix        \\
				        -ph=$temp_buffer_ph_map        \\
				        -n=$job_counter                \\
				        -gwas=$INPUT_GWAS_BLOCKS";
				        
	open (my $fh,  q{>}, $sge_script) or die("Unable to open $sge_script : $!");
	print $fh "#!/bin/bash\n";
	print $fh $command, "\n";
	close $fh;				        
	
	system "nice -5 qsub -v $qsub_opt_v -e $qsub_err -o $qsub_out -q newnodes.q $sge_script";
	#unlink $sge_script;	
	$job_counter++;	
	
	#$temp_counter--;
	#last if($temp_counter == 0);	        		
}
#check all jobs are finished
my $COMPLETE;
my $JOB_STATUS = `qstat`;
print $JOB_STATUS . "\n";

if ($JOB_STATUS =~ /\Q$sge_job_string/){
	$COMPLETE=0;
}else{
	$COMPLETE=1;
}
while($COMPLETE == 0){
	print "Waiting 10 seconds for $sge_job_string jobs to complete..\n";
	sleep 10;
	$JOB_STATUS = `qstat`;
	if ($JOB_STATUS =~ /\Q$sge_job_string/){
		$COMPLETE=0;
	}else{
		$COMPLETE=1;
	}			
}
print "All $sge_job_string jobs completed.\n";

#$qsub_out should contain the tsv data for this phenotype 
#eg
#2 hour glucose	1	0.079	1.85356811862836	0.89030513517407	0.0779220779220779
################
#7 gather all results strings in a unique data structure
################

#you need to open and close all files of the form VET_*.out
my %ph_to_stats; #key phenotype - value stats

my $file_string = $sge_job_string . '*.out';
my @files = <VET_JOBS_*.out>;
#my @files = <$file_string>;

unless ( (scalar @files) eq (scalar keys %phenotype_map)  ){
	my $a = scalar @files;
	my $b = scalar keys %phenotype_map;
	print STDERR "Attention: number of out files ($a) and phenotypes tested ($b) differ\n";
}

chdir $directory;
foreach my $file (@files) {
	my $data = `cat $file`;
	chomp $data;
	my ($item, $obs, $exp, $fc, $l2fc, $pval) = split("\t", $data);
	$ph_to_stats{$item} =  join("\t", $obs, $exp, $fc, $l2fc, $pval);
}
system "rm $QSUB_BASENAME*.*"; #rm qsub temps

#get pvalues to correct
my @pvals;
foreach my $item (sort keys %ph_to_stats){
        my @fields = split("\t", $ph_to_stats{$item});
        push(@pvals, $fields[4]);
}

##################
#8 call R and adjust pvalue:
###################

#I think p.adjust is more suitable than "qvalue" see posts here:
#https://stat.ethz.ch/pipermail/bioconductor/2014-January/056982.html
my $R = Statistics::R->new();
$R->set('pvals', \@pvals);
$R->run(q`adjustp <- p.adjust(pvals, method = "BH", n = length(pvals))`);
my $corrected_pvals = $R->get('adjustp');
$R->stop();

my $counter = 0;
open (my $outstream,  q{>}, $out_results) or die("Unable to open $out_results : $!");
print $outstream "TRAIT\tOBSERVED\tEXPECTED\tFOLD\tLOG2_FOLD\tPVAL\tADJ_PVAL\n";
foreach my $item (sort keys %ph_to_stats){
	print $outstream $item . "\t" . $ph_to_stats{$item} . "\t" . $$corrected_pvals[$counter] . "\n";
	$counter++;
}
close $outstream;
print "*FINISHED*\n";
unlink $out_temp_buffer;
unlink $out_temp_buffer_ph_list;
unlink $out_temp_bg_matrix;
unlink $temp_buffer_ph_map;








#------------------------------------------------------------------------------
#subs
#------------------------------------------------------------------------------
sub get_DAF{
	my ($hash) = @_;
	
	tie *FILE,   'IO::Zlib', $INPUT_BACKGROUND, "rb";
	while (<FILE>)	{ 
		chomp;
		next if($_ eq '');
		next if($_ =~ /^#/);
		#my $ALTF_AFR; #not doing anything with the african frequencies atm 
		my $AF_EUR; #ALTERNATE allele frequencies; many rare background SNPs don't have EUR allele freq
		my $AF_GLOB;
		my $ANC;
		my $FLAG; # set to one if the ancestral is the alternate	
	
		my @fields = split("\t", $_);
		#skip indels
		next if(length($fields[3]) > 1);
		next if(length($fields[4]) > 1);
		my $key =  $fields[0] . '-' . $fields[1];
		
		#get derived allele frequencies
		my $ref = $fields[3]; 
		my $alt = $fields[4];
		my @info = split(";", $fields[7]);
		
		#get ancestral allele info and allele frequencies=====================
		foreach my $item (@info){
			if($item =~ /^AA=(.*)/){
				$ANC = $1;
			}
			if($item =~ /^EUR_AF=(.+)/){
				$AF_EUR = $1;	
			}
			if($item =~ /^AF=(.+)/){
				$AF_GLOB = $1;
			}	
		}	
		next unless($ANC);
		
		if(uc($ANC) eq uc($alt)){
			$FLAG = 1;		
		}elsif(uc($ANC) eq uc($ref)){			
		}elsif( ($ANC eq '') or ($ANC eq 'N') or ($ANC eq '.') or ($ANC eq '-') ){
			next;
		}else{ 
			next;			
		}
	
		if($AF_GLOB){
			if($FLAG){
				if($AF_GLOB == 1){
						$$hash{$key}{$AFGLOB} = '0.0';			
				}else{
					$$hash{$key}{$AFGLOB} = (1 - $AF_GLOB);				
				}
			}else{
				$$hash{$key}{$AFGLOB} = $AF_GLOB;
			}
		}
		if($AF_EUR){
			if($FLAG){
				if($AF_EUR == 1){
						$$hash{$key}{$AFEUR} = '0.0';			
				}else{
					$$hash{$key}{$AFEUR} = (1 - $AF_EUR);				
				}
			}else{
				$$hash{$key}{$AFEUR} = $AF_EUR;
			}
		}		
	}
	close FILE;	
	return;
}


sub get_MAF{
	my ($hash) = @_;
	
	tie *FILE,   'IO::Zlib', $INPUT_BACKGROUND, "rb";
	while (<FILE>)	{ 
		chomp;
		next if($_ eq '');
		next if($_ =~ /^#/);
		#my $ALTF_AFR; #not doing anything with the african frequencies atm 
		my $AF_EUR; #ALTERNATE allele frequencies; many rare background SNPs don't have EUR allele freq
		my $AF_GLOB;
	
		my @fields = split("\t", $_);
		#skip indels
		next if(length($fields[3]) > 1);
		next if(length($fields[4]) > 1);
		my $key =  $fields[0] . '-' . $fields[1];
		
		my @info = split(";", $fields[7]);
		
		#get ancestral allele info and allele frequencies=====================
		foreach my $item (@info){
			if($item =~ /^EUR_AF=(.+)/){
				$AF_EUR = $1;	
			}
			if($item =~ /^AF=(.+)/){
				$AF_GLOB = $1;
			}	
		}	
		
		if($AF_GLOB){
			$$hash{$key}{$AFGLOB} = $AF_GLOB;
		}
		if($AF_EUR){
			$$hash{$key}{$AFEUR} = $AF_EUR;
		}		
	}
	close FILE;	
	return;
}

############### binning ############
sub bin_variants{
	my ($file, $hash, $type) = @_;

	my $skipped = 0;
	my $valid = 0;
	my $total = 0;
	my %local_hash; #for the background the vcf is not unique, so fill hash. Do it also for the foreground, it won't hurt

	tie *FILE,   'IO::Zlib', $file, "rb";

	while (<FILE>)	{
		chomp;
		next if($_ eq '');
		next if($_ =~ /^#/);
		my $chr; my $pos;
	
		if($type eq 'fg'){ #fg is a bed
			($chr, $pos) = (split /\t/)[0,2];
		}elsif($type eq 'bg'){
			($chr, $pos) = (split /\t/)[0,1]; #vcf file
		}else{
			print STDERR "bin_variants(): type: $type not recognised. Aborting\n..";
			exit -1;
		}	
		my $key = $chr . '-' . $pos;
		$local_hash{$key} = 1;
	}
	close FILE;

	foreach my $item (keys %local_hash){
		$total++;		
		#TODO DISCARDING those VDR-BVs which don't have clear ANCESTRAl/DERIVED info - see what you're wasting here
		unless($variants_1kg{$item}){
			#print STDERR "Warning: background VDR-BV variant $item: key not found in background. Skipping..\n";
			$skipped++;
			next;
		}
	#	unless($variants_1kg{$item}{$AFEUR}){
	#		print STDERR "Warning: background VDR-BV variant at coordinates $item does not have DAF CEU info. Skipping..\n";
	#		$skipped_bg_vars++;
	#		next;		
	#	}
		unless($variants_1kg{$item}{$AFGLOB}){
			#print STDERR "Warning: background VDR-BV variant at coordinates $item does not have GLOB MAF info. Skipping..\n";
			$skipped++;
			next;		
		}
		#test that the DAF value is numeric?
		my $bin_id = int( $variants_1kg{$item}{$AFGLOB} * $BIN_NUMBER ) / $BIN_NUMBER;
		#Here we care about the coordinate of the variants composing the bin, and each should be present once, because we suppose they're unique
		if($type eq 'fg'){
			$$hash{$bin_id}++;
		}else{
			$$hash{$bin_id}{$item}++; 
		}
		$valid++;
	}
	return ($total, $skipped, $valid);
}


###############histogram##############
sub print_histogram{
	my ($hash, $type) = @_;
	
	if($type eq 'fg'){
		for (0..$BIN_NUMBER){
	    	$_ = $_/$BIN_NUMBER;
	    	if($$hash{$_}){ 
			    my $number_in_this_bin = $$hash{$_};    	
		    	print STDERR "$_\t" .$number_in_this_bin . "\n";    		
	    	}else{
		    	print STDERR "$_\t0\n";	    		
	    	}
	    }		
	}elsif($type eq 'bg'){
		for (0..$BIN_NUMBER){
	    	$_ = $_/$BIN_NUMBER;
	    	if(ref $$hash{$_}){ 
     			my $number_in_this_bin = scalar(keys %{$$hash{$_}});    	
	    		print STDERR "$_\t" .$number_in_this_bin . "\n";   		
	    	}else{
		    	print STDERR "$_\t0\n";	    		
	    	}
	    }				
	}else{
	    	print STDERR "Error: print_histogram() - type $type not recognised. Aborting..\n";
	    	exit -1;		
	}
	return;	
}
