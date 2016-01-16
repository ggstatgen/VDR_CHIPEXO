#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#get all the rsID for the positions in the "interestinghet_consensus.txt" output of do_alleleseq_pp_01_intersect_results.pl
#this is like the other one, but I won't manipulate the pvalues, because the consensus file does not have them (each entry has at least 1 pvalue)


#map them from the original .vcf you used. Maybe grep is the fastest program?
#you will need these to send them to SNAP/Haploreg, get the snps in LD, check for disease association with the do_SNAP_... scripts.
#optionally, only select HET snps in motifs

#Since I'm testing collated variants, I will have to use both variant files

my $INPUT_VARIANTS_BASE = '/net/isi-scratch/giuseppe/VDR/VARIANTS/1000g_hg19_b37/ALL.phase1_release_v3.20101123.snps_indels_svs.genotypes_FIFTEEN_1K_GEN.vcf.gz';
my $INPUT_VARIANTS_IMPUTE = '/net/isi-scratch/giuseppe/VDR/VARIANTS/IMPUTATION/d_IMPUTE_out_chksize_4mb/IMPUTE_1kg_to_hapmap_autosomes_g1k_v37.vcf.gz';
#print STDERR "Using variant set in $INPUT_VARIANTS\n\n";


#structure:
###fileformat=VCFv4.1
###FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype calls">
###FORMAT=<ID=GP,Number=3,Type=Float,Description="Genotype call probabilities">
##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA06986 NA06989 NA06997 NA07029 NA07045 NA10831 NA10846 NA10847 NA11829 NA11832 NA11918 NA11919 NA12264 NA12383 NA12489 NA12716 NA12752 NA12872 NA19189 NA19190 NA19191
#chr1    10583   rs58108140      G       A       .       .       .       GT:GP   0/1:0.031,0.963,0.006   0/1:0.032,0.961,0.006   ./.:0.733,0.226,0.041   ./.:0.784,0.206,0.01    ./.:0.468,0.393,0.139   ./.:0.755,0.24,0.005    ./.:0.841,0.1

#skip the headers (##, #) and get chr pos id

#OUTPUT
#new file interestinghets (containing only associations in peaks) and containing one line with rsID(s)
#comma separated list of rsIDs
#I SKIP OVERLAPPING VARIANTS (within the sample, same position, multiple rsIDs)
#because these were skipped when building the personalised genome
#I count them and output the count to stderr.

#OPTIONAL INPUT:
#bed with weeder VDR interval for strong motifs
#you could only produce output for the intersection

my $infile_alleleseq;
my $infile_motifs; # file obtained by converting pscanchip ris to bed using "do_pscanchip_out_to_bed.pl"
my $outfile;
my $outfile_list;#comma separated list of rsIDs, useful as a direct input to HAPLOREG
my $allsites; #flag. If set, gets ASB sites even when they don't fall in a peak bed interval
my $qval;
GetOptions(
        'i=s'        =>\$infile_alleleseq,
        'allsites'   =>\$allsites,
        'qval=f'     =>\$qval,
        'm=s'        =>\$infile_motifs
);

#$infile_alleleseq = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ_15_SINGLES/PIPELINE_SANDBOX/interestingHets_NA06989.txt";
#$sample_id = "NA06989";

if(!$infile_alleleseq){
	print "USAGE: do_alleleseq_get_rsids_from_interestinghets.pl -i=<INFILE_ASSOC> -allsites -qval=<QVAL> -m=<INFILE_MOTIFS>\n";
	print "<INFILE_ASSOC> file interestingHets_overlaps.txt from do_alleleseq_pp_02..\n";
	print "(opt)<allsites> flag; if set, retrieve asb sites even when they did not fall in peak interval file.\n";
	print "(opt)<QVAL> qval threshold (if you want to restrict the returned sites, eg 0.05)\n";
	print "(opt)<INFILE_MOTIFS> file bed containing motif instances. Only HET associations intersecting these will be kept.\n";
    exit 1;
}

print STDERR "LIMITING OUTPUT TO MOTIF INTERVALS: $infile_motifs..\n" if($infile_motifs);

my ($filename, $directory) = fileparse($infile_alleleseq);
$filename =~ s/(.*)\..*/$1/;
#if($infile_motifs){
#	$outfile = $directory . "\/". $filename . '_snpids_bds_motifs.txt';
#}else{
#	$outfile = $directory . "\/". $filename . '_snpids_bds.txt';
#}
#$outfile_list = $directory . "\/". $filename . '_snpids_rsIDlist.txt';
my %rsIDlist;

#fill hash with motifs if any
my %motif_data;
if($infile_motifs){
	open (my $instream,  q{<}, $infile_motifs) or die("Unable to open $infile_motifs : $!");
	while(<$instream>){
		chomp;
		next if($_ =~ /^#/); 
		next if($_ =~ /^track type/);
		my ($chr, $start, $stop, $motif, $score) = (split /\t/)[0,1,2,3,4];
		my $coord = $start . '-' .  $stop;
		$motif_data{$chr}{$coord}{'motif'} = $motif;
		$motif_data{$chr}{$coord}{'score'} = $score;
	}
	close $instream;
}

#map the actual instances, remove associations not in binding sites (unless allsites), output enriched file
my $header;
my %output_data;
my $SKIPPED = 0;
open (my $instream,  q{<}, $infile_alleleseq) or die("Unable to open $infile_alleleseq : $!");
while(<$instream>){
	my $M_OVERLAP;
	my $vcf_line;
	
	chomp;
	next if($_ eq '');
	if($_ =~ /^chrm/){
		$header = $_;
		next;
	}
	my ($a_chr,$a_snppos,$ref,$a_bindingSite) = (split /\t/)[0,1,2,3];
	next if(!$a_chr);
	next if(!$a_snppos);
	
	#binding sites yes/no?
	unless($allsites){
		next unless ($a_bindingSite =~ /1/);
	}
	#$a_chr =~ s/(.+)/chr$1/; #hg19
	
	if($infile_motifs){
		#check for motif overlap
		my $c_chr = 'chr' . $a_chr;
		foreach my $coord_pair (sort keys %{ $motif_data{$c_chr} }){
			my ($m_start,$m_stop) = split("-",$coord_pair);
			if( ($a_snppos >= $m_start) && ($a_snppos <= $m_stop) ){
				$M_OVERLAP = 1;
				$output_data{$a_chr}{$a_snppos}{'motif'} = $motif_data{$c_chr}{$coord_pair}{'motif'};
				$output_data{$a_chr}{$a_snppos}{'score'} = $motif_data{$c_chr}{$coord_pair}{'score'};
				last;
			}
		}		
		next unless($M_OVERLAP);
	}
	
	#test chr and location against 1kg catalog, first the normal, then the impute if nothing found
	
	my $lines = `zgrep -P \"$a_chr\t$a_snppos\" $INPUT_VARIANTS_BASE`;
	if(!$lines){
		$lines = `zgrep -P \"$a_chr\t$a_snppos\" $INPUT_VARIANTS_IMPUTE`;
	}
	
	#lines will have to be split and $rsid captured
	my @lines = split("\n", $lines);
	if(scalar @lines eq 1){
		$vcf_line = $lines[0];
	}elsif(scalar @lines > 1){
		#IF we are here, there are two variants at the same position.
		#in my experience, this is when an indel overlaps a snp, where the indel is on one allele and the snp on the other
		#When vcf2diploid builds the personalised genomes, it SKIPS these. Therefore these variants are not incorporated 
		#in the personalised genomes
		#I skip them here as well, but count them.
		print STDERR "WARNING: position $a_chr - $a_snppos corresponds to more than one VCF line. Skipping..\n";
		$SKIPPED += 1;
		next;
	}else{
		print STDERR "WARNING: position $a_chr - $a_snppos not found in any of the two 1kg vcfs. Skipping..\n";
		next;
	}
	my @vcf_fields = split("\t", $vcf_line);
	my $snp_id = $vcf_fields[2];
	
	if(  (!$snp_id) or ($snp_id eq '')  ){
		print  STDERR "WARNING: I did not find an ID for this: $a_chr - $a_snppos. Skipping..\n";
		next;
	}
	$output_data{$a_chr}{$a_snppos}{'alleleseq_data'} =  $snp_id . "\t" . $_;
	$rsIDlist{$snp_id} = 1;	
}
close $instream;

print STDERR "Skipped overlapping positions: " . $SKIPPED . "\n";

#print output, sorted by chr followed by position
#open (my $outstream,  q{>}, $outfile) or die("Unable to open $outfile : $!");
if($infile_motifs){
	print STDOUT 'SNP_ID' . "\t" . $header . "\t" . 'MOTIF' . "\t" . 'MOTIF_SCORE' . "\n";
	foreach my $chr (sort keys %output_data){
		foreach my $snp_pos (sort keys %{ $output_data{$chr} }){
			print STDOUT $output_data{$chr}{$snp_pos}{'alleleseq_data'} . "\t" . 
				  $output_data{$chr}{$snp_pos}{'motif'}   . "\t" . 
				  $output_data{$chr}{$snp_pos}{'score'}   . "\n";
		}
	}
}else{
	print STDOUT 'SNP_ID' . "\t" . $header . "\n";
	foreach my $chr (sort keys %output_data){
		foreach my $snp_pos (sort keys %{ $output_data{$chr} }){
			print STDOUT $output_data{$chr}{$snp_pos}{'alleleseq_data'} . "\n";
		}
	}
}
#close $outstream;

my @rsIDs = keys %rsIDlist;
#open ($outstream,  q{>}, $outfile_list) or die("Unable to open $outfile_list : $!");
print STDERR join(",", @rsIDs) . "\n";
#close $outstream;
