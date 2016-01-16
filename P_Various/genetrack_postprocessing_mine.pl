#!/usr/bin/perl
use strict;
use warnings;
#use Bio::Tools::GFF;
use Getopt::Long;


#genetrack postprocessing uses bedtools merge
#bedtools merge merges features regardless of strand
#I want to ONLY merge features if they're at most 35bp apart and on different strand

#GFF format
#Fields are: <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
#from here http://www.sanger.ac.uk/resources/software/gff/spec.html

#INPUT: gff raw out of genetrack
#OUTPUT: gff without singletons and with paired peaks


#SAMPLE INPUT
##gff-version 3
#chr1    genetrack       .       564480  564490  30.3116900285   +       .       readcount=28;ID=P564485;stddev=3.03382294023;height=30.3116900285
#chr1    genetrack       .       564555  564565  34.104478709    +       .       readcount=37;ID=P564560;stddev=3.27640701853;height=34.104478709
#chr1    genetrack       .       564589  564599  48.9043737066   +       .       readcount=40;ID=P564594;stddev=2.51197133742;height=48.9043737066
#chr1    genetrack       .       564640  564650  38.3398207915   +       .       readcount=35;ID=P564645;stddev=1.97039310281;height=38.3398207915
#chr1    genetrack       .       564664  564674  43.4731681704   +       .       readcount=39;ID=P564669;stddev=2.8508854179;height=43.4731681704
#chr1    genetrack       .       564675  564685  45.7547328299   +       .       readcount=39;ID=P564680;stddev=2.56922565126;height=45.7547328299
#chr1    genetrack       .       564703  564713  28.7019951503   +       .       readcount=28;ID=P564708;stddev=2.72460444611;height=28.7019951503
#chr1    genetrack       .       564736  564746  39.3745670075   +       .       readcount=32;ID=P564741;stddev=2.7613402543;height=39.3745670075
#..
#chr1    genetrack       .       10031   10041   12.8811920398   -       .       readcount=13;ID=P9996;stddev=0.498518515262;height=12.8811920398
#chr1    genetrack       .       564502  564512  20.6386278892   -       .       readcount=21;ID=P564467;stddev=3.17158586769;height=20.6386278892
#chr1    genetrack       .       564525  564535  39.2902256588   -       .       readcount=41;ID=P564490;stddev=3.51033148388;height=39.2902256588
#chr1    genetrack       .       564584  564594  31.1662365808   -       .       readcount=32;ID=P564549;stddev=3.47634304083;height=31.1662365808
#chr1    genetrack       .       564607  564617  29.6117466725   -       .       readcount=24;ID=P564572;stddev=2.21382963813;height=29.6117466725
#chr1    genetrack       .       564620  564630  33.8855529973   -       .       readcount=34;ID=P564585;stddev=2.67323657233;height=33.8855529973
#chr1    genetrack       .       564660  564670  41.9819966876   -       .       readcount=38;ID=P564625;stddev=2.78101996076;height=41.9819966876
#chr1    genetrack       .       564684  564694  49.2238582827   -       .       readcount=43;ID=P564649;stddev=2.16616391111;height=49.2238582827
#chr1    genetrack       .       564713  564723  48.9210413545   -       .       readcount=53;ID=P564678;stddev=3.2219432764;height=48.9210413545
#chr1    genetrack       .       564747  564757  50.8566878988   -       .       readcount=53;ID=P564712;stddev=2.94569334649;height=50.8566878988
#chr1    genetrack       .       564770  564780  32.4298062429   -       .       readcount=31;ID=P564735;stddev=2.80181300008;height=32.4298062429


#remove sigletons

my $infile;
my $tagdist; # tag distance for merging gff signals
my $BEDTOOLS = '/net/isi-scratch/giuseppe/tools/bedtools-2.17.0/bin';

GetOptions(
	'i=s'  =>\$infile,
	'd=s' =>\$tagdist,
);
if(!$infile){
     print "USAGE: genetrack_postprocessing.pl -i=<INFILE> -d=<TAGDIST>\n";
     print "<INFILE> input .gff file\n";
     print "<TAGDIST> max tag distance to use to consider two +/- peaks as part of the same (e.g. 35)\n";
     exit 1;
}
if(!$tagdist){
     print "USAGE: genetrack_postprocessing.pl -i=<INFILE> -d=<TAGDIST>\n";
     print "<INFILE> input .gff file\n";
     print "<TAGDIST>  max tag distance to use to consider two +/- peaks as part of the same (e.g. 35)\n";
     exit 1;
}
my $file_nosingletons = $infile;
$file_nosingletons =~ s/(.*)\..*/$1/;

#my $file_nosingletons_plus = $file_nosingletons . "_nosgtons_plus.gff";
#my $file_nosingletons_minus = $file_nosingletons . "_nosgtons_minus.gff";
#my $file_nosingletons_plus_s = $file_nosingletons . "_nosgtons_plus_s.gff";
#my $file_nosingletons_minus_s = $file_nosingletons . "_nosgtons_minus_s.gff";
my $file_nosingletons_sorted = $file_nosingletons . "_nosgtons_s.gff";
my $file_nosingletons_paired = $file_nosingletons . "_nosgtons_p.bed";
my $file_final_bed = $file_nosingletons . "_nosgtons_final.bed";
#my $file_nosingletons_merged = $file_nosingletons . ".bed";
$file_nosingletons = $file_nosingletons . "_nosgtons.gff";

open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
open (my $outstream,  q{>}, $file_nosingletons) or die("Unable to open $file_nosingletons : $!");
while(<$instream>){
	#remove 0.0, 0.1, 0.2, 0.3 peaks
	print $outstream $_ unless($_ =~ /stddev=0\.[0|1|2|3]/);
	#print $outstream $_ unless($_ =~ /stddev=0\.0;/);
	#if($_ =~ /stddev=0\.[0|1|2|3]/){
#		print $_;
#	}else{
#		print $outstream $_;
#	}
}
close $instream;
close $outstream;

#sort -k1,1V -k2,2g 
system "sort -k1,1V -k4n,4 $file_nosingletons > $file_nosingletons_sorted";
open ($instream,  q{<}, $file_nosingletons_sorted) or die("Unable to open $file_nosingletons_sorted : $!");
open ($outstream, q{>}, $file_nosingletons_paired) or die("Unable to open $file_nosingletons_paired : $!");


#assume data is ordered
#at any time, keep two memory one buffers: one for + and one for minus
#for each line, get strand
#look at buffer 
#I merge only one + with one -. Therefore I only do the merging once I meet one type of strand. For example +
#If peak is
#1 unpaired
#2 has low number of reads
#3 low sd
#----------------remove - most likely duplication 
my $bufferline_plus;
my $bufferline_minus;
while(<$instream>){
	#print $outstream $_ if($_ =~ /^\#\#/);
	next if($_ =~ /^\#\#/);
	my $bed_line;
	my @strand_plus_fields;
	my @strand_minus_fields;
	my $c_minus_start;
	my $c_minus_end;
	my $c_plus_start;
	my $c_plus_end;
	my $chr;
	
	my $strand = (split /\t/)[6];
	#----------
	#the current data line is +
	#----------	
	if($strand eq '+'){
		if(!$bufferline_minus){
			#this is the first line
			$bufferline_plus = $_;
			next;
		}
		@strand_plus_fields = split(/\t/, $_);
		$chr =  $strand_plus_fields[0];
		$c_plus_start = $strand_plus_fields[3];
		$c_plus_end = $strand_plus_fields[4];
		#my $c_plus_start = (split /\t/)[3];
		#my $c_plus_end   = (split /\t/)[4];
		#check coordinates of last minus strand. Was it close enough for this to be considered a paired peak?
		@strand_minus_fields = split(/\t/, $bufferline_minus);
	    $c_minus_start = $strand_minus_fields[3];
		$c_minus_end = $strand_minus_fields[4];
		#are this + line and the former - line a pair?
		#If so, print the merged pair
		if( abs($c_plus_start - $c_minus_end) <= $tagdist ){
			#merge and create bed entry line
			$c_minus_start -= 1;
			$bed_line = "$chr\t$c_minus_start\t$c_plus_end\tm";
			print $outstream $bed_line, "\n";
			$bufferline_plus = $_;
			next;
		}else{
			#if this line is not the second in a pair, two things can happen:
			#1 line is the first in a peak
			#2 line is a singleton
			#in both cases, you can only tell it at the next iteration
			#so buffer it 
			#also, you need to deal with the former line as well
			$c_minus_start -= 1;
			$bed_line = "$chr\t$c_minus_start\t$c_minus_end\t-";
			print $outstream $bed_line, "\n";
			$bufferline_plus = $_;
			next;
		}
		
	#----------
	#the current data line is -
	#----------	
	}elsif($strand eq '-'){
		if(!$bufferline_plus){#this is the first line
			$bufferline_minus = $_;
			next;
		}
		@strand_minus_fields = split(/\t/, $_);
		$chr =  $strand_minus_fields[0];
		$c_minus_start = $strand_minus_fields[3];
		$c_minus_end = $strand_minus_fields[4];		
		#check coordinates of last plus strand. Was it close enough for this to be considered a paired peak?
		@strand_plus_fields = split(/\t/, $bufferline_plus);
		$c_plus_start = $strand_plus_fields[3];
		$c_plus_end = $strand_plus_fields[4];	
		
		#are this - line and the former + line a pair?
		#If so, print the merged pair
		if( abs($c_minus_start - $c_plus_end) <= $tagdist ){
			#merge and create bed entry line 
			$c_plus_start -= 1;
			$bed_line = "$chr\t$c_plus_start\t$c_minus_end\tm";
			print $outstream $bed_line, "\n";
			$bufferline_minus = $_;
			next;
		}else{
			#if this line is not the second in a pair, two things can happen:
			#1 line is the first in a peak
			#2 line is a singleton
			#in both cases, you can only tell it at the next iteration
			#so buffer it 
			#do some checks to remove single
			$bufferline_minus = $_;
			$c_plus_start -= 1;
			$bed_line = "$chr\t$c_plus_start\t$c_plus_end\t-";
			print $outstream $bed_line, "\n";
			next;
		}	
	}else{
		print "Unable to interpret strand info for line: $_\n";
		print "Aborting now\n";
		exit -1;		
	}
}

close $instream;
close $outstream;
#sort -k1,1V -k2,2g
#system "sort -k1,1V -k2,2g $file_nosingletons_paired > $file_nosingletons_paired_sorted";
#system "$BEDTOOLS/mergeBed -n -d $tagdist -i $file_nosingletons_sorted > $file_nosingletons_merged";
#system "$BEDTOOLS/mergeBed -n -i $file_nosingletons_paired_sorted > $file_final_bed";
system "sort -k1,1V -k2,2g $file_nosingletons_paired | $BEDTOOLS/bedtools merge -n -i stdin > $file_final_bed";
unlink $file_nosingletons_sorted;
unlink $file_nosingletons_paired;




#create file for plus strand and file for minus strand
#open ($instream,  q{<}, $file_nosingletons) or die("Unable to open $file_nosingletons : $!");
#open (my $out_plus, q{>},  $file_nosingletons_plus ) or die("Unable to open $file_nosingletons_plus  : $!");
#open (my $out_minus, q{>}, $file_nosingletons_minus) or die("Unable to open $file_nosingletons_minus : $!");
#while(<$instream>){
#	if($_ =~ /^\#\#/){
#		print $out_plus $_;
#		print $out_minus $_;
#		next;
#	}
#	
#	my $strand = (split /\t/)[6];
#	if($strand eq '+'){
#		print $out_plus $_;
#	}elsif($strand eq '-'){
#		print $out_minus $_;
#	}else{
#		print "Unable to assign strand: field six contains $_\n";
#		exit -1;
#	}
#}
#close $instream;
#close $out_plus;
#close $out_minus;

#system "sort -k1,1 -k4n,4 $file_nosingletons_plus > $file_nosingletons_plus_s";
#system "sort -k1,1 -k4n,4 $file_nosingletons_minus > $file_nosingletons_minus_s";
#unlink $file_nosingletons_plus;
#unlink $file_nosingletons_minus;
#
#open the two strand specific files
#open (my $instream_plus,  q{<}, $file_nosingletons_plus_s) or die("Unable to open $file_nosingletons_plus_s : $!");
#open (my $instream_minus,  q{<}, $file_nosingletons_minus_s) or die("Unable to open $file_nosingletons_minus_s : $!");
#
#my $data_plus  = read_file_line($instream_plus);
#my $data_minus = read_file_line($instream_minus);


###################
#sub read_file_line{
#  my $fh = shift;
#
#  if ($fh and my $line = <$fh>) {
#    chomp $line;
#    return [ split(/\t/, $line) ];
#  }
#  return;
#}


