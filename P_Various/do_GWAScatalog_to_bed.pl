#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use String::Approx 'amatch';

my $TOOL_PATH       = '/net/isi-scratch/giuseppe/tools';
my $BEDTOOLS  = $TOOL_PATH . '/bedtools2-2.19.0/bin';
my $CHROM_SIZE_PATH = '/net/isi-scratch/giuseppe/indexes/chrominfo/hg19.chrom.sizes';

#INFO
#THIS uses all annotation in GWAScatalog as is provided the term appears for at least x SNPs.

#INFO 20/11
#This converts the GWAS catalog .txt table found here
#http://www.genome.gov/gwastudies/ (The SNP data in the catalog have been mapped to dbSNP Build 132 and Genome Assembly, GRCh37/hg19. )
#in a bed file to use for GAT annotation

#SAMPLE INPUT
#Date Added to Catalog	PUBMEDID	First Author	Date	Journal	Link	Study	Disease/Trait	Initial Sample Size	Replication Sample Size	Region	Chr_id	Chr_pos	Reported Gene(s)	Mapped_gene	Upstream_gene_id	Downstream_gene_id	Snp_gene_ids	Upstream_gene_distance	Downstream_gene_distance	Strongest SNP-Risk Allele	SNPs	Merged	Snp_id_current	Context	Intergenic	Risk Allele Frequency	p-Value	Pvalue_mlog	p-Value (text)	OR or beta	95% CI (text)	Platform [SNPs passing QC]	CNV
#11/13/2013	23696099	Ding K	05/20/2013	G3 (Bethesda)	http://www.ncbi.nlm.nih.gov/pubmed/23696099	Genetic variants that confer resistance to malaria are associated with red blood cell traits in African-Americans: an electronic medical record-based genome-wide association study.	Red blood cell traits	1904 Afican American individuals	411 Afican American individuals	11p15.4	11	5230907	NA	OR51V1 - HBB	283111	3043		        8.98	       15.79	rs7120391-C	rs7120391	0	7120391	Intergenic	1	0.12	5E-9	8.301029995663981	(MCHC)	    .30	[0.20-0.40] unit increase	Illumina [907,954]	N

#FIELDs
#0 Date Added to Catalog
#1 PUBMEDID
#2 First Author
#3 Date
#4 Journal
#5 Link
#6 Study
#7 Disease/Trait
#8 Initial Sample Size
#9 Replication Sample Size
#10 Region
#11 Chr_id
#12 Chr_pos
#13 Reported Gene(s)
#14 Mapped_gene
#15 Upstream_gene_id
#16 Downstream_gene_id
#17 Snp_gene_ids
#18 Upstream_gene_distance
#19 Downstream_gene_distance
#20 Strongest SNP-Risk Allele
#21 SNPs
#22 Merged
#23 Snp_id_current
#24 Context
#25 Intergenic
#26 Risk Allele Frequency
#27 p-Value
#28 Pvalue_mlog
#29 p-Value (text)
#30 OR or beta
#31 95% CI (text)
#32 Platform [SNPs passing QC]
#33 CNV

#VALUE
#0 11/13/2013
#1 23696099
#2 Ding K
#3 05/20/2013
#4 G3 (Bethesda)
#5 http://www.ncbi.nlm.nih.gov/pubmed/23696099
#6 Genetic variants that confer resistance to malaria are associated with red blood cell traits in African-Americans: an electronic medical record-based genome-wide association study.
#7 Red blood cell traits
#8 1904 Afican American individuals
#9 411 Afican American individuals
#10 11p15.4
#11 11
#12 5230907
#13 NA
#14 OR51V1 - HBB
#15 283111
#16 3043	
#17
#18 8.98
#19 15.79
#20 rs7120391-C
#21 rs7120391
#22 0
#23 7120391
#24 Intergenic
#25 1
#26 0.12
#27 5E-9
#28 8.301029995663981
#29 (MCHC)
#30 .30
#31 [0.20-0.40] unit increase
#32 Illumina [907,954]
#33 N

#I will need the description field (7), chromosome, position

#ATM I'll pick the following field numbers:
#7 11	12 21


my $infile;
my $min_snp_number;
my $SLOP;
GetOptions(
        'i=s'      =>\$infile,
        'minsnp=i' =>\$min_snp_number,
        'slop=i'   =>\$SLOP
);
if(!$infile){
     print "USAGE: do_GWAScatalog_to_bed.pl -i=<INFILE> -minsnp=<MIN_SNP_NUMBER> -slop=<BASEPAIRS>\n";
     print "<INFILE> GWAScatalog txt file\n";
     print "<MIN_SNP_NUMBER> minimum number of SNP per phenotype to consider for output\n";
     print "<BASEPAIRS> number of bp requested upstream (downstream) of each SNP (default: 150,000)\n";
     print "MAKE SURE the gwascatalog is mapped to hg19 before continuing\n";
     exit 1;
}
if(!$min_snp_number){
     print "USAGE: do_GWAScatalog_to_bed.pl -i=<INFILE> -minsnp=<MIN_SNP_NUMBER> -slop=<BASEPAIRS>\n";
     print "<INFILE> GWAScatalog txt file\n";
     print "<MIN_SNP_NUMBER>  minimum number of SNP per phenotype to consider for output\n";
     print "<BASEPAIRS> number of bp requested upstream (downstream) of each SNP\n";
     exit 1;
}
$SLOP = '150000' if(!$SLOP);

my ($filename, $directory) = fileparse($infile);
$filename =~ s/(.*)\..*/$1/;
my $outfile_segments = $directory . "\/" .  $filename . '_segments.bed';

#1st pass
#collect unique entries (study,disease/trait,chromosome,position,snp id) in hash
open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
my %unique_entries;
while(<$instream>){
	chomp;
	next if($_ eq '');
	next if($_ =~ /^Date/); #header

	#There are duplicate SNP-phenotype lines.
	my ($study, $dis_trait, $chrom, $pos, $snp_id) = (split /\t/)[6,7,11,12,21];
	next if(!$dis_trait);
	next if($dis_trait eq '');
	next if($snp_id eq '');
	next if(!$chrom);
	next if($chrom eq '');
	next if(!$pos);
	next if($pos eq '');
	
	my $line = $study . "\t" . $dis_trait . "\t" . $chrom . "\t" . $pos . "\t" . $snp_id;
#	print $line, "\n";
#exit;

	if(!$unique_entries{$line}){
		$unique_entries{$line} = 1;
	}else{
		$unique_entries{$line} += 1;		
	}
}
close $instream;


#2nd pass
#build phenotype->snp-number map
my %phenotype_map;
foreach my $entry (keys %unique_entries){
	my ($study, $dis_trait, $chrom, $pos, $snp_id) = split (/\t/, $entry);
	if(!$phenotype_map{lc($dis_trait)}){
		$phenotype_map{lc($dis_trait)} = 1;
	}else{
		$phenotype_map{lc($dis_trait)} += 1;		
	}
}



#create one bed for each of the terms above, if the term has > $min_snp_number
my %collective_bed;
foreach my $item (sort keys %phenotype_map){
	next if($phenotype_map{$item} < $min_snp_number);
	
	my $outfile_single_phenotype           = $directory . "\/" .  $filename . '_phenotype_' . '_temp.bed';
	my $outfile_single_phenotype_processed = $directory . "\/" .  $filename . '_phenotype_' . '_temp_out.bed';
	open (my $temp_stream,  q{>}, $outfile_single_phenotype) or die("Unable to open $outfile_single_phenotype : $!");
	foreach my $entry (keys %unique_entries){
		my ($study, $dis_trait, $chrom, $pos, $snp_id) = split (/\t/, $entry);
		if(lc($dis_trait) eq $item){
			if($chrom =~ /23/){
				$chrom = 'chrX';
			}elsif($chrom =~ /24/){
				$chrom = 'chrY';
			}elsif($chrom =~ /25/){
				$chrom = 'chrM';
			}elsif($chrom =~ /\d+/){
				$chrom = 'chr' . $chrom;
			}else{
				print "Chromosome string $chrom is not recognised. Skipping entry..\n";
				next;
			}
			#use 0-based start coordinates, which are the proper BED format (pos-1, pos)
			my $line = $chrom . "\t" . ($pos-1)  . "\t" . $pos . "\t" . $item . "\n";
			print $temp_stream 	$line;
		}
	}
	close $temp_stream;
	#you have a bed with only the snp and its  location. Extend location to +/- 150.000, sort, merge
	system "$BEDTOOLS/bedtools slop -b $SLOP -i $outfile_single_phenotype -g $CHROM_SIZE_PATH |  sort -k1,1V -k2,2g | $BEDTOOLS/bedtools merge -i stdin > $outfile_single_phenotype_processed";
	open (my $temp_stream_processed,  q{<}, $outfile_single_phenotype_processed) or die("Unable to open $outfile_single_phenotype_processed : $!");
	while(<$temp_stream_processed>){
		chomp;
		next if($_ eq '');
		my $data_line = $_ . "\t" . $item . "\n";
		$collective_bed{$data_line} = 1;
	}
	close $temp_stream_processed;
	unlink $outfile_single_phenotype;
	unlink $outfile_single_phenotype_processed;
}
open (my $segment_stream,  q{>}, $outfile_segments) or die("Unable to open $outfile_segments : $!");
foreach my $item (sort keys %collective_bed){ print $segment_stream $item; }
close $segment_stream;
