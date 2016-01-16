#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#snp to peak associations
#I obtain a matrix of TMM normalised counts from the bam files using Diffbind (EdgeR)
#This matrix is one of the imputs for SNPTEST and the scripts do_association_SNPTEST_*.pl
#However it needs some postprocessing

#1 add a # at the beginning of the header
#2 rename header samples GM > NA
#3 remove three columns: 
#4 insert columns of 0s for missing samples

#4 is important: the matrix must contain one column FOR EACH OF THE THIRTY SAMPLES.
#The reason is due to how I imputed my samples and the .bgen format of IMPUTE2. It DOES NOT LABEL the sample columns (1 0 0 0 1 0 ecc) so it WON'T KNOW IF A SAMPLE IS MISSING
#THEREFORE create fake columns if you decide to avoid using some samples (and then create an exclusion list for those in snptest)

#inputs
#1 ordered sample list (hardcoded)
#2 matrix out of diffbind

my $infile_rcmatrix;
GetOptions(
        'm=s'        =>\$infile_rcmatrix
);

#$infile_rcmatrix = "/net/isi-scratch/giuseppe/VDR/POOLED_4/08_BAM_STAMPY_g1k_v37_MIN21/d_macs2_nodup_q0.05/d_bed_RBL/d_diffbind/vdr_o3_matrix.txt";

if(!$infile_rcmatrix){
     print "USAGE: do_association_diffbind_postprocess_matrix.pl  -m=<RC_MATRIX>\n";
     print "<RC_MATRIX> raw matrix of read counts produced by DiffBind. row: peak. column: sample\n";
     exit 1;
}
my($basename, $directory) = fileparse($infile_rcmatrix);
$basename =~ s/(.*)\..*/$1/;
my $outfile_rcmatrix = $directory . $basename . '_p.txt';
my $outfile_exclusion_file = $directory . $basename . '_exclusion.txt';

my %FULL_SAMPLE_LIST = ("NA06986" => 1, 
				   "NA06989" => 1,
				   "NA06997" => 1,
				   "NA07029" => 1,
				   "NA07045" => 1, 
				   "NA10831" => 1, 
				   "NA10846" => 1, 
				   "NA10847" => 1, 
				   "NA11829" => 1, 
				   "NA11832" => 1, 
				   "NA11918" => 1, 
				   "NA11919" => 1, 
				   "NA12264" => 1, 
				   "NA12383" => 1, 
				   "NA12489" => 1, 
				   "NA12716" => 1, 
				   "NA12752" => 1, 
				   "NA12872" => 1, 
				   "NA19189" => 1, 
				   "NA19190" => 1, 
				   "NA19191" => 1,
				   "NA19213" => 1,
				   "NA19214" => 1,
				   "NA19215" => 1,
				   "NA19235" => 1,
				   "NA19236" => 1,
				   "NA19237" => 1,
				   "NA19247" => 1,
				   "NA19248" => 1,
				   "NA19249" => 1);

my $HEADER ="#\tNA06986\tNA06989\tNA06997\tNA07029\tNA07045\tNA10831\tNA10846\tNA10847\tNA11829\tNA11832\tNA11918\tNA11919\tNA12264\tNA12383\tNA12489\tNA12716\tNA12752\tNA12872\tNA19189\tNA19190\tNA19191\tNA19213\tNA19214\tNA19215\tNA19235\tNA19236\tNA19237\tNA19247\tNA19248\tNA19249\n";

#sample 
#CHR     START   END     GM06986 GM06989 GM06997 GM07029 GM07045 GM10831 GM10846 GM10847 GM11829 GM11832 GM11918 GM11919 GM12264 GM12383 GM12489 GM12716 GM12752 GM12872 GM19189 GM19190 GM19191 GM19213 GM19214 GM19215 GM19235 GM19236 GM19237 GM19247 GM19248 GM19249
#1       1       9998    10035   5.68865708820131        4.93210402291583        3.05874603054557        2.2892937304082 5.79892942893578        3.32679811586275        #6.55865646637907        4.73704907404362        1.93295146585962
#        6.56610716955441        3.28234297220948        4.02812009303887        5.06848647729742        2.63187642682116        2.8088436085977 3.55801161513053        #10.7798361576846        1.15871376452785        1.43020113902865
#        2.07466244821304        3.62963732132452        2.34395715205964        3.11894479732402        8.24870484363639        7.85403979553074        3.262079664881  #6.42764461654137        1.90477123273068        2.66263151553064
#        0.99294940912575

#todo 
#get samples from matrix header
#check if anyone is missing based on hash
#if sample in missing, put a NA or 0 in corresponding column
#compare header to full header
#replace header
#remove column 1,2,3


my @sample_names;
my %sample_names;
open (my $instream,      q{<}, $infile_rcmatrix) or die("Unable to open $infile_rcmatrix : $!");
open (my $outstream,     q{>}, $outfile_rcmatrix) or die("Unable to open $outfile_rcmatrix : $!");
while(<$instream>){
	chomp;
	my %sample_name_to_rd;
	#header
	#get available samples, print new header, return
	if($_ =~ /^\CHR/){
		#replace GM with NA in sample names
		(my $header = $_) =~ s/GM/NA/g;
		@sample_names = split("\t", $header);
		@sample_names = @sample_names[ 3 .. $#sample_names ];
		foreach my $item (@sample_names){ $sample_names{$item} = 1; }
		
		#check which files are missing from this set and fill exclusion file
		open (my $outstream_ex,     q{>}, $outfile_exclusion_file) or die("Unable to open $outfile_exclusion_file : $!");
		foreach my $item (sort keys %FULL_SAMPLE_LIST){ print $outstream_ex $item, "\n" unless($sample_names{$item}); }
		close $outstream_ex;
		#print header
		print $outstream $HEADER;
		next;
	}
	
	#get a structure with sample name to count
	my @data = split("\t", $_);
	my $peak_id = shift @data;
	my @counts = @data[3 .. $#data];
	@sample_name_to_rd{@sample_names} = @counts;
	
	#now you print the data
	#No matter how many columns there are, you print 31 columns: peak_n rd_sample1 ... rd_sample30
	
	print $outstream $peak_id . "\t";
	foreach my $item (sort keys %FULL_SAMPLE_LIST){
		if($sample_name_to_rd{$item}){
			print $outstream $sample_name_to_rd{$item}, "\t";
		}else{
			#print $outstream 'NA', "\t"; #compatible with SNPTEST, but max @counts then fails..
			print $outstream '0.0', "\t";
			
		}
	}
	print $outstream "\n";
}
close $instream;
close $outstream;