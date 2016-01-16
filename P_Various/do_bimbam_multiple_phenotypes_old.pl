#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#SLOW: here I did one intersection per iteration. Superseded by do_bimbam_multiple_phenotypes.pl

#14/3/2014
#carry out bayesian regression in a snp-peak association analysis setting 
#things to try:
#1 different peakset
#2 different vcf set
#3 different bimbam model?
#4 multiple snps around one peak? (say +/- 200bp?)

#bimbam for bayes factor based SNP-phenotype association
#https://www.bcm.edu/research/centers/childrens-nutrition-research-center/mcmcmc/index.cfm?pmid=18981

#at the moment this file assumes you have generated the read counts and normalised them using DESeq2.

#INPUTS: 
#1 vcf file with snps and genotype
#2 bed file with peak set
#3 matrix of counts for that peakset

#process:
#for each peak, intersect peak with vcf and get all the snps intersecting. Create BIMBAM INPUTs and run bimbam. Get snps-position to BF in hash
#EXTREMELY slow. Alternative? Can I do only one bedtools intersection and examine it one line at a time?

#current pipeline:
#1 get the matrix of counts in a data structure indexed by peak and sample
#2 intersect vcf and peaks once with bedtools
#3 using the output of 2, build datastructure to keep all the data:
#peak
# ---chr
# ---start
# ---end
# ---variation
# ------snpID
# ---------ref
# ---------alt
# ---------gt
# ---------readcount?


my $infile_vcf;
my $infile_bed;
my $infile_rcmatrix;
#eg /net/isi-scratch/giuseppe/VDR/POOLED_4/07_BAM_STAMPY_MIN21/d_htseq_counts/d_gold_26/normalized_counts.txt
GetOptions(
        'vcf=s'         =>\$infile_vcf,
        'bed=s'         =>\$infile_bed,
        'rc=s'          =>\$infile_rcmatrix
);
if(!$infile_vcf){
     print "USAGE: do_bimbam_multiple_phenotypes.pl -vcf=<VCF_FILE> -bed=<BED_FILE> -rc=<READ_COUNT>\n";
     print "<VCF_FILE> initial set of variation_genotypes to use (eg hapmap, omni2.5, 1000g)\n";
     print "<BED_FILE> set of peaks to use for the association analysis\n";
     print "<READ_COUNT> matrix of normalised read counts for the peaks in <BED_FILE> + headers. row: peak. column: sample\n";
     exit 1;
}
if(!$infile_bed){
     print "USAGE: do_bimbam_multiple_phenotypes.pl -vcf=<VCF_FILE> -bed=<BED_FILE> -rc=<READ_COUNT>\n";
     print "<VCF_FILE> initial set of variation_genotypes to use (eg hapmap, omni2.5, 1000g)\n";
     print "<BED_FILE> set of peaks to use for the association analysis\n";
     print "<READ_COUNT> matrix of normalised read counts for the peaks in <BED_FILE> + headers. row: peak. column: sample\n";
     exit 1;
}
if(!$infile_rcmatrix){
     print "USAGE: do_bimbam_multiple_phenotypes.pl -vcf=<VCF_FILE> -bed=<BED_FILE> -rc=<READ_COUNT>\n";
     print "<VCF_FILE> initial set of variation_genotypes to use(eg hapmap, omni2.5, 1000g)\n";
     print "<BED_FILE> set of peaks to use for the association analysis\n";
     print "<READ_COUNT> matrix of normalised read counts for the peaks in <BED_FILE> + headers. row: peak. column: sample\n";
     exit 1;
}
my $TOOL_PATH       = '/net/isi-scratch/giuseppe/tools';
my $BEDTOOLS  = $TOOL_PATH . '/bedtools2-2.19.0/bin';
my $BIMBAM    = $TOOL_PATH . '/bimbam';
my($basename, $directory) = fileparse($infile_rcmatrix);
$basename =~ s/(.*)\..*/$1/;

my $global_out = $directory . '/' . $basename . 'bimbam_out.txt';

#TEMP FILES
my $infile_bed_intersect_vcf    = $directory . 'TEMP_bimbam_bed_intersect_vcf.txt';

my $temp_peak_file              = $directory . 'TEMP_bimbam_regression_single.bed';
my $temp_vcf_bed_intersect_file = $directory . 'TEMP_bimbam_regression_single.vcf';
my $temp_phenotype_file         = $directory . 'TEMP_bimbam_regression_ph.txt';
my $temp_genotype_file          = $directory . 'TEMP_bimbam_regression_gt.txt';
my $temp_snploc_file            = $directory . 'TEMP_bimbam_regression_snploc.txt';

#############
#1. get the matrix of read counts in a data structure indexed by peak and sample
# also get a hash of candidate sample names (a subset of those in the vcf for which you want to computer regression)
#############
my @sample_names; 
my %peak_to_counts;
my %candidate_sample_names; #samples for which we have read counts, could be a subset of those in the vcf
open (my $instream,     q{<}, $infile_rcmatrix) or die("Unable to open $infile_rcmatrix : $!");
while(<$instream>){
	chomp;
	my %sample_name_to_rd;
	if($_ =~ /^\#/){
		@sample_names = split("\t", $_);
		shift @sample_names;
		foreach (@sample_names){ s/\"//g; }
		%candidate_sample_names = map { $_ => 1 } @sample_names;
		next;
	}
	my @counts = split("\t", $_);
	my $this_peak_ID = shift @counts;
	$this_peak_ID =~ s/\"//g;
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
my @cols;
open ($instream,     q{<}, $infile_vcf) or die("Unable to open $infile_vcf : $!");
while(<$instream>){
	chomp;
	next if ($_ =~ /^\#\#/); #metadata	
	if($_ =~ /^\#CHROM/){ #header
			@cols = split("\t", $_);
		    @c_names = @cols[9..(@cols-1)];
		    next;
	}
}
close $instream;


###############################################OLD#########################################

 

########
#3. cycle bed, peak by peak
########
my %bimbam_data_structure;
#in the hash above you'll save the collective data and print it in a file
#you could make it similar to the bimbam out "single"
#rsnum           pos     chr     bf      se      rank    pv      mu      a       d
#rs1539637       1080925 1       -0.013   NA      1               4.24e-01               1.724   -0.059  0.000

open ($instream,     q{<}, $infile_bed) or die("Unable to open $infile_bed : $!");
while(<$instream>){ #main loop
	next if($_ =~ /^\#/); #header
	open (my $temp_bed_stream,     q{>}, $temp_peak_file) or die("Unable to open $temp_peak_file : $!");
	#get name - field 3
	my $peak_name = (split /\t/)[3];
	print $temp_bed_stream $_;
	close $temp_bed_stream;
	
	#you have a peak and its ID
	########
	#4. get intersection between vcf and bed, if no snp intersects this peak, move on
	########
	system "$BEDTOOLS/bedtools intersect -a $infile_vcf -b $temp_peak_file > $temp_vcf_bed_intersect_file";
	if(-z $temp_vcf_bed_intersect_file){
		unlink $temp_bed_stream;
		unlink $temp_vcf_bed_intersect_file;
		next;		
	} 
	
	################
	#5. create phenotype file from the reads data structure
	################
	open (my $temp_phenotype_stream,     q{>}, $temp_phenotype_file) or die("Unable to open $temp_phenotype_file : $!");
	foreach my $sample (sort keys %{ $peak_to_counts{$peak_name} }){
		print $temp_phenotype_stream  $peak_to_counts{$peak_name}{$sample}, "\n";
	}
	close $temp_phenotype_stream;
	#foreach my $rc (@reads){ print $temp_phenotype_stream $rc, "\n"; }

	
	#We have a peak with at least a snp
	################
	#6. Collect snp(s), genotypes and coords in data structure
	#	this won't contain more than 2/3 snps
	################
	my $snp_counter = `cat $temp_vcf_bed_intersect_file | wc -l`;
	chomp $snp_counter;
	my %sample_to_gt;
	my %snp_collection;
	open (my $temp_vcf_stream,     q{<}, $temp_vcf_bed_intersect_file) or die("Unable to open $temp_vcf_bed_intersect_file : $!");
	#what to do in case of multiple snps?
	while(<$temp_vcf_stream>){
		chomp;
		#bedtools should put no metadata in the intersection, but check
		next if ($_ =~ /^\#\#/); #metadata
		next if ($_ =~ /^\#CHROM/); #metadata
		
		my @snp_line = split(/\t/, $_);
		#collect coordinates
		my $snp_chr = $snp_line[0];
		next if(!$snp_chr);
		my $snp_pos = $snp_line[1];
		next if(!$snp_pos);
		my $snp_id  = $snp_line[2];
		next if(!$snp_id);
		next if ($snp_id eq '.');
		my $snp_ref  = $snp_line[3];
		next if(!$snp_ref);
		next if($snp_ref =~ /^[ACGT]{2,}/);
		my $snp_alt  = $snp_line[4];
		next if(!$snp_alt);
		next if($snp_alt =~ /^[ACGT]{2,}/);
		
		#unused at the moment
		#my $snp_info = $snp_line[7];
		#my $snp_ancestral = get_ancestral_allele($snp_info); # NOT USED see get_sample_to_genotype_hash
		
		#collect genotypes
		my @genotypes = @snp_line[9..(@cols-1)]; 
		@sample_to_gt{@c_names} = @genotypes; # for some of these I won't have normalised reads
		#the next function assumes there is ONLY ONE ALTERNATE ALLELE (not a list)
		#if this is not the case, next
		if($snp_alt =~ /(\w\W)+/){
			print "ERROR: alternate allele string: $snp_alt seem to contain multiple alternate alleles. Skipping..\n";
			next;
		}
		my $s_to_g = get_sample_to_genotype_hash($snp_ref, $snp_alt, %sample_to_gt);
		$snp_collection{$snp_id}{'chr'} =  $snp_chr;
		$snp_collection{$snp_id}{'pos'} =  $snp_pos;
		foreach my $sample (sort keys % {$s_to_g} ) {
			$snp_collection{$snp_id}{'gt'}{$sample} = $$s_to_g{$sample};
		}	
	}
	close $temp_vcf_stream;
	
	################
	#7. create bimbam genotype file, genotype location file
	################
	open (my $temp_gt_stream,     q{>}, $temp_genotype_file) or die("Unable to open $temp_genotype_file : $!");
	print $temp_gt_stream @sample_names . "\n"; #this needs to be obtained from the rd file
	print $temp_gt_stream $snp_counter, "\n";
	print $temp_gt_stream 'IND, ' . join(", ", @sample_names) . "\n";
	foreach my $this_snp_id ( sort keys %snp_collection ){
		print $temp_gt_stream $this_snp_id;
		foreach my $this_sample (@sample_names){
			print $temp_gt_stream ', ' . $snp_collection{$this_snp_id}{'gt'}{$this_sample};
		}
#		foreach my $this_sample ( sort keys %{ $snp_collection{$this_snp_id}{'gt'} } ){
#			print $temp_gt_stream ', ' . $snp_collection{$this_snp_id}{'gt'}{$this_sample};
#		}
		print $temp_gt_stream "\n";
	}
	close $temp_gt_stream;
	
	################
	#8. create bimbam genotype location file
	################
	open (my $temp_snploc_stream,  q{>}, $temp_snploc_file) or die("Unable to open $temp_snploc_file : $!");
		foreach my $this_snp_id ( sort keys %snp_collection ){
			#turn chromosome into number
			my $chr = $snp_collection{$this_snp_id}{'chr'};
			if($chr =~ /^chr(\d+)/){
				$chr = $1;
			}elsif($chr =~ /^chrX/i){
				$chr = '23';
			}elsif($chr =~ /^chrY/i){
				$chr = '24';
			}elsif($chr =~ /^chrM/i){
				$chr = '25';
			}else{
				print "Unrecognised chromosome string: $chr. Aborting.\n";
				exit;
			}
			print $temp_snploc_stream $this_snp_id . ', ' . $snp_collection{$this_snp_id}{'pos'} . ', ' . $chr . "\n";
		}
	close $temp_snploc_stream;	
	
	################
	#10. run bimbam
	################	
	#typical options
	#-g bimbam_hapmapr27_genotypes_hg19_INTERSECT_vdr_o3_CONSENSUS_gold26_gt.txt
	#-p vdr_o3_consensus_normalizedcounts_07bamstampy_gold26.txt
	#-pos bimbam_hapmapr27_genotypes_hg19_INTERSECT_vdr_o3_CONSENSUS_snploc.txt
	#-sort (sorts the output snp by BF)
	#-pval 1000
	#-o test1
	my $out_prefix = 'bimbam_peak_' . $peak_name;
	system "$BIMBAM/bimbam_lin -g $temp_genotype_file -p $temp_phenotype_file -pos $temp_snploc_file -pval 1000  -sort -o $out_prefix";
	
	################
	#11. open outdata and collect BF and pvalue. Keep in datastructure
	################
	#the data will be in /output/$out_prefix.single.txt
	#format:
	#\#\# bf=log10_BF  se=log10_std_err        pv=p-value      
	#rsnum           pos     chr     bf      se      rank    pv      mu      a       d
	#rs1539637       1080925 1       -0.013   NA      1               4.24e-01               1.724   -0.059  0.000
	my $bimbam_outfile         = $directory . '/output/' .  $out_prefix . '.single.txt';
	open (my $bimbam_output_stream,     q{<}, $bimbam_outfile) or die("Unable to open $bimbam_outfile : $!");
	while(<$bimbam_output_stream>){
		chomp;
		next if ($_ =~ /^\#\#/);
		next if ($_ =~ /^rsnum/);
		my ($snp_id, $snp_pos, $chr, $logBF, $logP, $mu, $a, $d) = (split /\t/)[0,1,2,3,6,7,8,9];
		$bimbam_data_structure{$chr}{$snp_pos} = join("\t", $snp_id, $logBF, $logP, $mu, $a, $d);
	}
	close $bimbam_output_stream;
	
	unlink $temp_peak_file;
	unlink $temp_vcf_bed_intersect_file;
	#bimbam inputs:
	unlink $temp_phenotype_file;
	unlink $temp_genotype_file;
	unlink $temp_snploc_file;
}
close $instream;

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

