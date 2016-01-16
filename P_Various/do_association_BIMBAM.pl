#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#14/3/2014
#carry out bayesian regression in a snp-peak association analysis setting 
#things to try:
#1 different peakset
#2 different vcf set (done - imputed)
#3 different bimbam model?
#4 multiple snps around one peak? (say +/- 200bp?)
#5 different normalisation for read count?

#bimbam for bayes factor based SNP-phenotype association
#https://www.bcm.edu/research/centers/childrens-nutrition-research-center/mcmcmc/index.cfm?pmid=18981

#at the moment this file assumes you have generated the read counts and normalised them using DiffBind.

#23/5/2014
#Implement quantile normalisation of the input vector using R

#INPUTS: 
#1 vcf file with snps and genotype
#2 bed file with peak set
#3 matrix of counts for that peakset

#current pipeline:
#1 get the matrix of counts in DATA_STRUCTURE_1 indexed by peak and sample
#2 intersect vcf and peaks ONCE with bedtools
#3 using the output of 2, build DATA_STRUCTURE_2 to keep all the data:
#peak
# ---chr
# ---start
# ---end
# ---variation
# ------snpID
# ---------ref
# ---------alt
# ---------gt
# ---------readcount? (this is still in DATA_STRUCTURE_1)

my $infile_vcf;
my $infile_bed;
my $infile_rcmatrix;
GetOptions(
        'vcf=s'         =>\$infile_vcf,
        'bed=s'         =>\$infile_bed,
        'rc=s'          =>\$infile_rcmatrix
);

#eg /net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_htseq_counts/d_gold_26/normalized_counts.txt
#$infile_vcf = '/net/isi-scratch/giuseppe/VDR/VARIANTS/omni25_hg19_b37/Omni25_genotypes_2141_samples.hg19_OMNI_25_GEN_PASS.vcf';
#$infile_vcf = '/net/isi-scratch/giuseppe/VDR/VARIANTS/1000g_hg19_b37/ALL.phase1_release_v3.20101123.snps_indels_svs.genotypes_FIFTEEN_1K_GEN_hg19.vcf';
#$infile_vcf = '/net/isi-scratch/giuseppe/VDR/VARIANTS/IMPUTATION/d_IMPUTE_out_chksize_4mb/IMPUTE_1kg_to_hapmap_autosomes_hg19.vcf';
#$infile_bed = '/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_diffbind/vdr_o3_peaks.bed';
#$infile_rcmatrix = '/net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21_RBL/d_CONSENSUS_PEAKSET_GEM_MACS/d_diffbind/vdr_o3_matrix_bimbam.txt';

if(!$infile_vcf){
     print "USAGE: do_ASSOCIATION_bimbam.pl -vcf=<VCF_FILE> -bed=<BED_FILE> -rc=<READ_COUNT>\n";
     print "<VCF_FILE> initial set of variation_genotypes to use (eg hapmap, omni2.5, 1000g)\n";
     print "<BED_FILE> set of peaks to use for the association analysis\n";
     print "<READ_COUNT> matrix of normalised read counts for the peaks in <BED_FILE> + headers. row: peak. column: sample\n";
     exit 1;
}
if(!$infile_bed){
     print "USAGE: do_ASSOCIATION_bimbam.pl -vcf=<VCF_FILE> -bed=<BED_FILE> -rc=<READ_COUNT>\n";
     print "<VCF_FILE> initial set of variation_genotypes to use (eg hapmap, omni2.5, 1000g)\n";
     print "<BED_FILE> set of peaks to use for the association analysis\n";
     print "<READ_COUNT> matrix of normalised read counts for the peaks in <BED_FILE> + headers. row: peak. column: sample\n";
     exit 1;
}
if(!$infile_rcmatrix){
     print "USAGE: do_ASSOCIATION_bimbam.pl -vcf=<VCF_FILE> -bed=<BED_FILE> -rc=<READ_COUNT>\n";
     print "<VCF_FILE> initial set of variation_genotypes to use(eg hapmap, omni2.5, 1000g)\n";
     print "<BED_FILE> set of peaks to use for the association analysis\n";
     print "<READ_COUNT> matrix of normalised read counts for the peaks in <BED_FILE> + headers. row: peak. column: sample\n";
     exit 1;
}
my $TOOL_PATH       = '/net/isi-scratch/giuseppe/tools';
my $BEDTOOLS  = $TOOL_PATH . '/bedtools2-2.20.1/bin';
my $BIMBAM    = $TOOL_PATH . '/bimbam';
my $RSCRIPT   = $TOOL_PATH . '/R-3.1.0/bin/Rscript';
my($basename, $directory) = fileparse($infile_rcmatrix);
$basename =~ s/(.*)\..*/$1/;

my $global_out = $directory . '/' . 'BIMBAM_' . $basename . '_global_out.txt';
#TEMP FILES
my $infile_bed_intersect_vcf    = $directory . 'TEMP_bimbam_bed_intersect_vcf.txt';
my $temp_phenotype_file         = $directory . 'TEMP_bimbam_regression_ph.txt';
my $temp_phenotype_file_qn      = $directory . 'TEMP_bimbam_regression_ph_qn.txt'; # quantile normalised phenotypes
my $temp_R_file                 = $directory . 'TEMP_bimbam_regression_ph.R';
my $temp_genotype_file          = $directory . 'TEMP_bimbam_regression_gt.txt';
my $temp_snploc_file            = $directory . 'TEMP_bimbam_regression_snploc.txt';

#############
#1. get the matrix of read counts in a data structure indexed by peak and sample
# also get a hash of candidate sample names (a subset of those in the vcf for which you want to computer regression)
#############
my %peak_to_counts;
my @sample_names;
my %candidate_sample_names; #samples for which we have read counts, could be a subset of those in the vcf
open (my $instream,     q{<}, $infile_rcmatrix) or die("Unable to open $infile_rcmatrix : $!");
while(<$instream>){
	chomp;
	my %sample_name_to_rd;
	if($_ =~ /^\#/){
		@sample_names = split("\t", $_);
		shift @sample_names;
		#foreach (@sample_names){ s/\"//g; } #DEseq matrix only
		%candidate_sample_names = map { $_ => 1 } @sample_names;
		next;
	}
	my @counts = split("\t", $_);
	my $this_peak_ID = shift @counts;
	#$this_peak_ID =~ s/\"//g; DEseq matrix only
	@sample_name_to_rd{@sample_names} = @counts;
	
	foreach my $s_name (keys %sample_name_to_rd){
		$peak_to_counts{$this_peak_ID}{$s_name} = $sample_name_to_rd{$s_name};
	}
}
close $instream;

#############
#2. get vcf metadata (name of samples for ALL the samples in the vcf)
#############
my @c_names;
open ($instream,     q{<}, $infile_vcf) or die("Unable to open $infile_vcf : $!");
while(<$instream>){
	chomp;
	my @cols;
	next if ($_ =~ /^\#\#/); #metadata	
	if($_ =~ /^\#CHROM/){ #header
			@cols = split("\t", $_);
		    @c_names = @cols[9..(@cols-1)];
		    last;
	}
}
close $instream;

########
#3. intersect peak file with vcf file
########
print "Intersecting snp file with bed file..\n";
#whatever the bed, only do the intersection on its first 4 fields: chr, start, stop, name
#this should help avoiding problems due to macs/gem/whatever adding different meta fields to the bed
system "cut -f 1,2,3,4 $infile_bed | $BEDTOOLS/bedtools intersect -wo -a stdin -b $infile_vcf > $infile_bed_intersect_vcf";
#output structure is:

##columns 0-3 (BED)
#chr1    1080816 1081095 6
##columns 4-12 (SNP)       
#chr1    1080925 rs1539637       G       C       .       .       AC=17;AN=20;SF=0,1      GT
##columns 13-? (GENOTYPES)     
#./.     ./.     ./.     1/1     ./.     1/1     1/1     1/1     1/1     1/1     ./.     ./.

###############
#4. create data structure
###############
my %peak_to_snp_to_genotype_map;
my %sample_to_gt;
my %final_sample_set;
open ($instream,     q{<}, $infile_bed_intersect_vcf) or die("Unable to open $infile_bed_intersect_vcf : $!");
while(<$instream>){
	my @peak_snp_line = split(/\t/, $_);
	pop(@peak_snp_line);#we don't need the last field produced by bedtools
	
	#collect peak data--------------------------------------------------------
	my $peak_chr   = $peak_snp_line[0];
	my $peak_start = $peak_snp_line[1];
	my $peak_end   = $peak_snp_line[2];
	my $peak_ID    = $peak_snp_line[3];
	unless($peak_to_snp_to_genotype_map{$peak_ID}){
		$peak_to_snp_to_genotype_map{$peak_ID}{'chr'} = $peak_chr;
		$peak_to_snp_to_genotype_map{$peak_ID}{'start'} = $peak_start;
		$peak_to_snp_to_genotype_map{$peak_ID}{'end'} = $peak_end;
	}
	
	#collect snp data---------------------------------------------------------
	my $snp_chr = $peak_snp_line[4];
	next if(!$snp_chr);
	if($snp_chr ne $peak_chr){
		print "WARNING: chromosome different for peak/snp pair: $snp_chr, $peak_chr. Skipping..\n";
		next;
	}
	my $snp_pos = $peak_snp_line[5];
	if( ($snp_pos < $peak_start) ||  ($snp_pos > $peak_end) ){
		print "WARNING: snp coordinates:$snp_pos out of peak coordinates boundaries: ($peak_start, $peak_end). Skipping..\n";
		next;
	}
	next if(!$snp_pos);
	my $snp_ID  = $peak_snp_line[6];
	next if(!$snp_ID);
	next if ($snp_ID eq '.');
	my $snp_ref  = $peak_snp_line[7];
	next if(!$snp_ref);
	next if($snp_ref =~ /^[ACGT]{2,}/); #we are excluding indels here
	my $snp_alt  = $peak_snp_line[8];
	next if(!$snp_alt);
	next if($snp_alt =~ /^[ACGT]{2,}/); #we are excluding indels here
	
	$peak_to_snp_to_genotype_map{$peak_ID}{'variation'}{$snp_ID}{'pos'} = $snp_pos;
	$peak_to_snp_to_genotype_map{$peak_ID}{'variation'}{$snp_ID}{'ref'} = $snp_ref;
	$peak_to_snp_to_genotype_map{$peak_ID}{'variation'}{$snp_ID}{'alt'} = $snp_alt;
	
	#collect genotype data------------------------------------------------------
	#does it collect indels?
	my @genotypes = @peak_snp_line[13..$#peak_snp_line]; 
	
	#map all samples in the vcf to their genotype:
	#(for some of these I won't have normalised reads)
	@sample_to_gt{@c_names} = @genotypes; 
	
	#the next function assumes there is ONLY ONE ALTERNATE ALLELE (not a list)
	#if this is not the case, next
	if($snp_alt =~ /(\w\W)+/){
		print "WARNING: alternate allele string: $snp_alt seem to contain multiple alternate alleles. Skipping..\n";
		next;
	}
	my $s_to_g = get_sample_to_genotype_hash($snp_ref, $snp_alt, %sample_to_gt);
	foreach my $sample (sort keys % {$s_to_g} ) {
		#only keep genotypes for samples in the read matrix file:
		next unless($candidate_sample_names{$sample});
		$final_sample_set{$sample} = 1 unless ($final_sample_set{$sample}) ;
		$peak_to_snp_to_genotype_map{$peak_ID}{'variation'}{$snp_ID}{'gt'}{$sample} = $$s_to_g{$sample};
	}
}
close $instream;

#clean up----------
undef %candidate_sample_names;
undef %sample_to_gt;
undef @c_names;
#------------------

#get array of final sample names (needed for genotype file)
my @final_sample_set;
foreach my $item (sort keys %final_sample_set) {
    push(@final_sample_set,$item);
}

################
#5 run genotype-phenotype association analysis
# Collect data from the two datastructures (peak-sample-reads, and peak-snp-genotype) and build temp bimbam files
################
my %bimbam_data_structure;
foreach my $peak_ID (sort keys %peak_to_snp_to_genotype_map){
	my $snp_counter = 0;
	#how many snps for this peak?
	$snp_counter++ foreach (sort keys %{ $peak_to_snp_to_genotype_map{$peak_ID}{'variation'}} );
	next if($snp_counter == 0);
	################
	#5.1 BIMBAM - quantile normalize phenotypes & create phenotype file 
	################
	open (my $temp_phenotype_stream,     q{>}, $temp_phenotype_file) or die("Unable to open $temp_phenotype_file : $!");
	foreach my $sample (sort keys %{ $peak_to_counts{$peak_ID} }){
		print $temp_phenotype_stream  $peak_to_counts{$peak_ID}{$sample}, "\n" if($final_sample_set{$sample});
	}
	close $temp_phenotype_stream;
	#quantile normalise using R
	open (my $temp_R_stream,     q{>}, $temp_R_file) or die("Unable to open $temp_R_file : $!");
	print $temp_R_stream 'x<-read.delim(\'' . $temp_phenotype_file . '\',header=F)' . "\n";
	print $temp_R_stream 'x <- x[,1]' . "\n";
	print $temp_R_stream 'x_qn = qqnorm(x,plot.it = F)$x' . "\n";
	print $temp_R_stream 'write.table(x_qn, file=\'' . $temp_phenotype_file_qn . '\', quote=F,row.names=F,col.names=F)' . "\n";
	close $temp_R_stream;
	system "$RSCRIPT $temp_R_file";
	#temp_phenotype_file_qn should contain the usable phenotypes
	
	################
	#5.2 BIMBAM - create genotype file
	################
	open (my $temp_gt_stream,     q{>}, $temp_genotype_file) or die("Unable to open $temp_genotype_file : $!");
	print $temp_gt_stream @final_sample_set . "\n"; #this needs to be obtained from the number of genotypes in the $peak_to_snp_to_genotype_map data structure
	print $temp_gt_stream $snp_counter  . "\n";
	print $temp_gt_stream 'IND, ' . join(", ", @final_sample_set) . "\n";
	foreach my $this_snp_id (sort keys %{ $peak_to_snp_to_genotype_map{$peak_ID}{'variation'}} ){
		print $temp_gt_stream $this_snp_id;
		foreach my $this_sample (sort keys %{ $peak_to_snp_to_genotype_map{$peak_ID}{'variation'}{$this_snp_id}{'gt'}} ){
			print $temp_gt_stream ', ' . $peak_to_snp_to_genotype_map{$peak_ID}{'variation'}{$this_snp_id}{'gt'}{$this_sample};
		}
		print $temp_gt_stream "\n";
	}
	close $temp_gt_stream;
	
	################
	#5.3 BIMBAM -  create genotype location file
	################
	open (my $temp_snploc_stream,  q{>}, $temp_snploc_file) or die("Unable to open $temp_snploc_file : $!");
		foreach my $this_snp_id (sort keys %{ $peak_to_snp_to_genotype_map{$peak_ID}{'variation'}} ){
			#turn chromosome into number
			my $chr = $peak_to_snp_to_genotype_map{$peak_ID}{'chr'};
			if( ($chr =~ /^chr(\d+)/) || ($chr =~ /^(\d+)/) ){
				$chr = $1;
			}elsif( ($chr =~ /^chrX/i) || ($chr =~ /^X/i) ){
				$chr = '23';
			}elsif( ($chr =~ /^chrY/i) || ($chr =~ /^Y/i) ){
				$chr = '24';
			}elsif( ($chr =~ /^chrM/i) || ($chr =~ /^M/i) ){
				$chr = '25';
			}else{
				print "Unrecognised chromosome string: $chr. Aborting.\n";
				exit;
			}
			print $temp_snploc_stream $this_snp_id . ', ' . $peak_to_snp_to_genotype_map{$peak_ID}{'variation'}{$this_snp_id}{'pos'} . ', ' . $chr . "\n";
		}
	close $temp_snploc_stream;
	
	################
	#5.4 run bimbam
	################	
	#typical options
	#-g bimbam_hapmapr27_genotypes_hg19_INTERSECT_vdr_o3_CONSENSUS_gold26_gt.txt
	#-p vdr_o3_consensus_normalizedcounts_07bamstampy_gold26.txt
	#-pos bimbam_hapmapr27_genotypes_hg19_INTERSECT_vdr_o3_CONSENSUS_snploc.txt
	#-sort (sorts the output snp by BF)
	#-pval 1000
	#-o test1
	
	my $out_prefix = 'bimbam_peak_' . $peak_ID;
	if($snp_counter >= 1){
		system "$BIMBAM/bimbam_lin -g $temp_genotype_file -p $temp_phenotype_file_qn -pos $temp_snploc_file -pval 10000 -o $out_prefix";
	}else{
		print "ERROR: no snps found in peak: $peak_ID. Skipping..\n";
		next;		
	}
	
	
	#multiple test: if there > 1 snps in the interval, run a multiple analysis 
#	if($snp_counter == 1){
#		system "$BIMBAM/bimbam_lin -g $temp_genotype_file -p $temp_phenotype_file_qn -pos $temp_snploc_file -pval 10000 -o $out_prefix";	
#
#	}elsif($snp_counter > 1 && $snp_counter <= 3){
#		system "$BIMBAM/bimbam_lin -g $temp_genotype_file -p $temp_phenotype_file_qn -pos $temp_snploc_file -pval 10000 -l $snp_counter -o $out_prefix";
#	}
#	elsif($snp_counter > 3){ #for speed reasons I test multi-snp bf for all subsets of size 1,2,3 only
#		system "$BIMBAM/bimbam_lin -g $temp_genotype_file -p $temp_phenotype_file_qn -pos $temp_snploc_file -pval 1000 -l 3 -o $out_prefix";
#	}
#	else{
#		print "ERROR: no snps found in peak: $peak_ID. Skipping..\n";
#		next;
#	}

	################
	#5.5 open outdata and collect BF and pvalue. Keep in datastructure
	################
	#the data will be in /output/$out_prefix.single.txt
	#format:
	#\#\# bf=log10_BF  se=log10_std_err        pv=p-value      
	#rsnum           pos     chr     bf      se      rank    pv      mu      a       d
	#rs1539637       1080925 1       -0.013   NA      1               4.24e-01               1.724   -0.059  0.000
	
	#this file will contain one row per snp.

	my $bimbam_outfile         = $directory . 'output/' .  $out_prefix . '.single.txt';
	open (my $bimbam_output_stream,     q{<}, $bimbam_outfile) or die("Unable to open $bimbam_outfile : $!");
	while(<$bimbam_output_stream>){
		chomp;
		next if ($_ =~ /^\#\#/);
		next if ($_ =~ /^rsnum/);
		#my $ciccia = $_;
		my ($bb_snp_id, $bb_snp_pos, $bb_chr, $bb_logBF, $bb_logP, $bb_mu, $bb_a, $bb_d) = (split /\t/)[0,1,2,3,7,9,10,11];
		next if($bb_chr eq '');
		next if($bb_snp_pos eq '');
		$bimbam_data_structure{$bb_chr}{$bb_snp_pos} = join("\t", $bb_snp_id, $bb_logBF, $bb_logP, $bb_mu, $bb_a, $bb_d);
	}
	close $bimbam_output_stream;
	
	#unlink temp files (bimbam inputs)
	unlink $temp_phenotype_file;
	unlink $temp_genotype_file;
	unlink $temp_snploc_file;
}


#print global output
#header
#rsnum           pos     chr     bf      se      rank    pv      mu      a       d
open (my $global_output_stream,     q{>}, $global_out) or die("Unable to open $global_out : $!");
print $global_output_stream "chr\tpos\trsnum\tbf\tpval\tmu\ta\td\n";
foreach my $chr (sort {$a<=>$b} keys %bimbam_data_structure){
	foreach my $pos ( sort {$a<=>$b} keys %{ $bimbam_data_structure{$chr} } ){
		print $global_output_stream $chr . "\t" . $pos . "\t" . $bimbam_data_structure{$chr}{$pos} . "\n";
	}
}
close $global_output_stream;

###########
#subs
###########

#the desired output is a hash with key the sample name, and value a string as follows:
#XX
#?? if no data available
#example input
#T C 1/1:0.8888      0/1:0.8888      1/1:0.8888      0/1:0.8888      0/1:0.8004 
sub get_sample_to_genotype_hash{
	my ($my_vcf_ref, $my_vcf_alt, %my_sample_to_gt) = @_;
	#my ($my_vcf_ref, $my_vcf_alt, $my_vcf_anc, %my_sample_to_gt) = @_;
	my %sample_to_genotype;

	my %allele_map = (
		'0' => $my_vcf_ref,
		'1' => $my_vcf_alt 
	);	

	foreach my $item (sort keys %my_sample_to_gt){
		my @genotype_fields = split(':', $my_sample_to_gt{$item});
		my $sample_gt = $genotype_fields[0];
		#check that the genotype is actually present: 
		if($sample_gt =~ /(\d{1})[\/|\|](\d{1})/){
			if(!$allele_map{$1}) { print "Error: genotype $sample_gt contains more than one alternate allele: $1. Aborting..\n"; }
			if(!$allele_map{$2}) { print "Error: genotype $sample_gt contains more than one alternate allele: $2. Aborting..\n"; }
			my $output = $allele_map{$1} . $allele_map{$2};
			$sample_to_genotype{$item} = $output;
		}else{
			#genotype no present
			$sample_to_genotype{$item} = 'NN';
		}
	}
	return \%sample_to_genotype;
}


################################################OLD#########################################
#
# 
#
#########
##3. cycle bed, peak by peak
#########
#my %bimbam_data_structure;
##in the hash above you'll save the collective data and print it in a file
##you could make it similar to the bimbam out "single"
##rsnum           pos     chr     bf      se      rank    pv      mu      a       d
##rs1539637       1080925 1       -0.013   NA      1               4.24e-01               1.724   -0.059  0.000
#
#open ($instream,     q{<}, $infile_bed) or die("Unable to open $infile_bed : $!");
#while(<$instream>){ #main loop
#	next if($_ =~ /^\#/); #header
#	open (my $temp_bed_stream,     q{>}, $temp_peak_file) or die("Unable to open $temp_peak_file : $!");
#	#get name - field 3
#	my $peak_name = (split /\t/)[3];
#	print $temp_bed_stream $_;
#	close $temp_bed_stream;
#	
#	#you have a peak and its ID
#	########
#	#4. get intersection between vcf and bed, if no snp intersects this peak, move on
#	########
#	system "$BEDTOOLS/bedtools intersect -a $infile_vcf -b $temp_peak_file > $temp_vcf_bed_intersect_file";
#	if(-z $temp_vcf_bed_intersect_file){
#		unlink $temp_bed_stream;
#		unlink $temp_vcf_bed_intersect_file;
#		next;		
#	} 
#	
#	################
#	#5. create phenotype file from the reads data structure
#	################
#	open (my $temp_phenotype_stream,     q{>}, $temp_phenotype_file) or die("Unable to open $temp_phenotype_file : $!");
#	foreach my $sample (sort keys %{ $peak_to_counts{$peak_name} }){
#		print $temp_phenotype_stream  $peak_to_counts{$peak_name}{$sample}, "\n";
#	}
#	close $temp_phenotype_stream;
#	#foreach my $rc (@reads){ print $temp_phenotype_stream $rc, "\n"; }
#
#	
#	#We have a peak with at least a snp
#	################
#	#6. Collect snp(s), genotypes and coords in data structure
#	#	this won't contain more than 2/3 snps
#	################
#	my $snp_counter = `cat $temp_vcf_bed_intersect_file | wc -l`;
#	chomp $snp_counter;
#	my %sample_to_gt;
#	my %snp_collection;
#	open (my $temp_vcf_stream,     q{<}, $temp_vcf_bed_intersect_file) or die("Unable to open $temp_vcf_bed_intersect_file : $!");
#	#what to do in case of multiple snps?
#	while(<$temp_vcf_stream>){
#		chomp;
#		#bedtools should put no metadata in the intersection, but check
#		next if ($_ =~ /^\#\#/); #metadata
#		next if ($_ =~ /^\#CHROM/); #metadata
#		
#		my @snp_line = split(/\t/, $_);
#		#collect coordinates
#		my $snp_chr = $snp_line[0];
#		next if(!$snp_chr);
#		my $snp_pos = $snp_line[1];
#		next if(!$snp_pos);
#		my $snp_id  = $snp_line[2];
#		next if(!$snp_id);
#		next if ($snp_id eq '.');
#		my $snp_ref  = $snp_line[3];
#		next if(!$snp_ref);
#		next if($snp_ref =~ /^[ACGT]{2,}/);
#		my $snp_alt  = $snp_line[4];
#		next if(!$snp_alt);
#		next if($snp_alt =~ /^[ACGT]{2,}/);
#		
#		#unused at the moment
#		#my $snp_info = $snp_line[7];
#		#my $snp_ancestral = get_ancestral_allele($snp_info); # NOT USED see get_sample_to_genotype_hash
#		
#		#collect genotypes
#		my @genotypes = @snp_line[9..(@cols-1)]; 
#		@sample_to_gt{@c_names} = @genotypes; # for some of these I won't have normalised reads
#		#the next function assumes there is ONLY ONE ALTERNATE ALLELE (not a list)
#		#if this is not the case, next
#		if($snp_alt =~ /(\w\W)+/){
#			print "ERROR: alternate allele string: $snp_alt seem to contain multiple alternate alleles. Skipping..\n";
#			next;
#		}
#		my $s_to_g = get_sample_to_genotype_hash($snp_ref, $snp_alt, %sample_to_gt);
#		$snp_collection{$snp_id}{'chr'} =  $snp_chr;
#		$snp_collection{$snp_id}{'pos'} =  $snp_pos;
#		foreach my $sample (sort keys % {$s_to_g} ) {
#			$snp_collection{$snp_id}{'gt'}{$sample} = $$s_to_g{$sample};
#		}	
#	}
#	close $temp_vcf_stream;
#	
#	################
#	#7. create bimbam genotype file, genotype location file
#	################
#	open (my $temp_gt_stream,     q{>}, $temp_genotype_file) or die("Unable to open $temp_genotype_file : $!");
#	print $temp_gt_stream @sample_names . "\n"; #this needs to be obtained from the rd file
#	print $temp_gt_stream $snp_counter, "\n";
#	print $temp_gt_stream 'IND, ' . join(", ", @sample_names) . "\n";
#	foreach my $this_snp_id ( sort keys %snp_collection ){
#		print $temp_gt_stream $this_snp_id;
#		foreach my $this_sample (@sample_names){
#			print $temp_gt_stream ', ' . $snp_collection{$this_snp_id}{'gt'}{$this_sample};
#		}
##		foreach my $this_sample ( sort keys %{ $snp_collection{$this_snp_id}{'gt'} } ){
##			print $temp_gt_stream ', ' . $snp_collection{$this_snp_id}{'gt'}{$this_sample};
##		}
#		print $temp_gt_stream "\n";
#	}
#	close $temp_gt_stream;
#	
#	################
#	#8. create bimbam genotype location file
#	################
#	open (my $temp_snploc_stream,  q{>}, $temp_snploc_file) or die("Unable to open $temp_snploc_file : $!");
#		foreach my $this_snp_id ( sort keys %snp_collection ){
#			#turn chromosome into number
#			my $chr = $snp_collection{$this_snp_id}{'chr'};
#			if($chr =~ /^chr(\d+)/){
#				$chr = $1;
#			}elsif($chr =~ /^chrX/i){
#				$chr = '23';
#			}elsif($chr =~ /^chrY/i){
#				$chr = '24';
#			}elsif($chr =~ /^chrM/i){
#				$chr = '25';
#			}else{
#				print "Unrecognised chromosome string: $chr. Aborting.\n";
#				exit;
#			}
#			print $temp_snploc_stream $this_snp_id . ', ' . $snp_collection{$this_snp_id}{'pos'} . ', ' . $chr . "\n";
#		}
#	close $temp_snploc_stream;	
#	
#	################
#	#10. run bimbam
#	################	
#	#typical options
#	#-g bimbam_hapmapr27_genotypes_hg19_INTERSECT_vdr_o3_CONSENSUS_gold26_gt.txt
#	#-p vdr_o3_consensus_normalizedcounts_07bamstampy_gold26.txt
#	#-pos bimbam_hapmapr27_genotypes_hg19_INTERSECT_vdr_o3_CONSENSUS_snploc.txt
#	#-sort (sorts the output snp by BF)
#	#-pval 1000
#	#-o test1
#	my $out_prefix = 'bimbam_peak_' . $peak_name;
#	system "$BIMBAM/bimbam_lin -g $temp_genotype_file -p $temp_phenotype_file -pos $temp_snploc_file -pval 1000  -sort -o $out_prefix";
#	
#	################
#	#11. open outdata and collect BF and pvalue. Keep in datastructure
#	################
#	#the data will be in /output/$out_prefix.single.txt
#	#format:
#	#\#\# bf=log10_BF  se=log10_std_err        pv=p-value      
#	#rsnum           pos     chr     bf      se      rank    pv      mu      a       d
#	#rs1539637       1080925 1       -0.013   NA      1               4.24e-01               1.724   -0.059  0.000
#	my $bimbam_outfile         = $directory . '/output/' .  $out_prefix . '.single.txt';
#	open (my $bimbam_output_stream,     q{<}, $bimbam_outfile) or die("Unable to open $bimbam_outfile : $!");
#	while(<$bimbam_output_stream>){
#		chomp;
#		next if ($_ =~ /^\#\#/);
#		next if ($_ =~ /^rsnum/);
#		my ($snp_id, $snp_pos, $chr, $logBF, $logP, $mu, $a, $d) = (split /\t/)[0,1,2,3,6,7,8,9];
#		$bimbam_data_structure{$chr}{$snp_pos} = join("\t", $snp_id, $logBF, $logP, $mu, $a, $d);
#	}
#	close $bimbam_output_stream;
#	
#	unlink $temp_peak_file;
#	unlink $temp_vcf_bed_intersect_file;
#	#bimbam inputs:
#	unlink $temp_phenotype_file;
#	unlink $temp_genotype_file;
#	unlink $temp_snploc_file;
#}
#close $instream;
#
##print global output
##header
##rsnum           pos     chr     bf      se      rank    pv      mu      a       d
#open (my $global_output_stream,     q{>}, $global_out) or die("Unable to open $global_out : $!");
#print $global_output_stream "chr\tpos\trsnum\tbf\tpval\tmu\ta\td\n";
#foreach my $chr (sort {$a<=>$b} keys %bimbam_data_structure){
#	foreach my $pos ( sort {$a<=>$b} keys %{ $bimbam_data_structure{$chr} } ){
#		print $global_output_stream $chr . "\t" . $pos . "\t" . $bimbam_data_structure{$chr}{$pos} . "\n";
#	}
#}
#close $global_output_stream;

#code to produce a manhattan plot with the results

