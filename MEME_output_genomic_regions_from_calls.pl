#!/usr/bin/perl
#Giuseppe  Oct 2013
#based on this script found here
#http://www-hsc.usc.edu/~valouev/QuEST/MEME_analysis_of_GABP.html
#That script was written by the QUEST writer and required a quest peak file which looks like this:

#R-<region_id> <chromosome> <region start>-<region end> ChIP: <ChIP enrichment intensity of the strongest peak> control: <control enrichment intensity at the position of the strongest ChIP peak> max_pos: <coordinate of the highest peak in the region> ef: <normalized enrichment fold at the position of the strongest peak> ChIP_tags: <ChIP reads starting within the region> background_tags: <background reads starting within the region> tag_ef: <normalized tag enrichment fold within the region> ps: <local peak shift estimate> cor: <local strand correlation value> os: <local oversampling value> z_score: <z-value for deviating from the Poisson expectation> n_lg_qv: <negative log_10 of Poisson p-value corrected for the number of tests>
#P-<region_id>-<peak_id> <chromosome> <peak coordinate> ChIP: <ChIP enrichment intensity> control: <control enrichment intensity> region: <region start>-<region end> ef: <normalized_enrichment fold> ChIP_tags: <ChIP reads starting with n the region> background_tags: <background reads starting within the region> tag_ef: <normalized tag enrichment fold within the region> ps: <local peak shift estimate> cor: <local strand correlation value> os: <local oversampling value> _score: NA n_lg_qv: <negative log_10 of approximated FDR>

#eg
#R-1 chr8 123640927-123641979 ChIP: 0.58 control: 0.0027 max_pos: 123641348 ef: 262.007 ChIP_tags: 660 background_tags: 16 tag_ef: 50.3 ps: 34 cor: 0.824 os: 1.54 z_score: 178.603 n_lg_qv: 40.731
#P-1-1 chr8 123641348 ChIP: 0.58 control: 0.0027 region: 123640927-123641979 ef: 262.007 ChIP_tags: 267 background_tags: 6 tag_ef: 54.3 ps: 34 cor: 0.832 os: 1.54 z_score: NA n_lg_qv: 2.118
#There are regions, these contain peaks, each peak belongs to a region

#eg MACS bed file
#chr1    1310782 1310939 MACS2_855D3_peak_1      92
#chr1    1590480 1590607 MACS2_855D3_peak_2      120
#chr1    2493424 2493611 MACS2_855D3_peak_3      663
#chr1    3582091 3582202 MACS2_855D3_peak_4      131
#chr1    3593518 3593645 MACS2_855D3_peak_5      94
#chr1    3713085 3713203 MACS2_855D3_peak_6      78
#eg GPS bed file
#chr1    894542  894742  1:894642           13.0 
#chr1    948852  949052  1:948952           42.0 
#chr1    968504  968704  1:968604           21.0 
#chr1    1057573 1057773 1:1057673          15.0 
#chr1    1167307 1167507 1:1167407          26.0 

#I want a version that:
#1 takes as input a generic bed file
#2 uses bedtools for the bed to fasta conversion 

use strict;
use warnings;
use diagnostics;

my $window_size = 100;

my $usage = qq{
    output_genomic_regions_from_calls.pl

    -----------------------------
    mandatory parameters:
    
    -i bed peak file
    -o output_file_prefix
    -r ref_genome file (should contain all contigs)

    -----------------------------
    optional parameters

    -w <window> default $window_size    

};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

## mandatory arguments
my $empty = "";
my $calls_fname = $empty;
my $output_fname_prefix = $empty;
my $ref_fname = $empty;
#my $output_fname_prefix = "";

## optional arguments
## parse command line arguments
while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-i') {$calls_fname = shift @ARGV;}
    elsif ( $this_arg eq '-o') {$output_fname_prefix = shift @ARGV;}
    elsif ( $this_arg eq '-r') {$ref_fname = shift @ARGV;}
    elsif ( $this_arg eq '-w') {$window_size = shift @ARGV;}

    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

print "\n";
print "Peaks        :           $calls_fname\n";
print "Ref genome   :           $ref_fname\n";
print "output prefix:           $output_fname_prefix\n";
print "window_size  :           $window_size\n";
print "\n";

if($calls_fname eq $empty || $output_fname_prefix eq $empty || $ref_fname eq $empty){
    print $usage;
}

open calls_file, "< $calls_fname" || die "calls_fname: $!\n";
my @peaks;
my $peaks_counter = 0;
while(<calls_file>){
	chomp;
    my $cur_line = $_;
	my ($cur_chrom, $region_begin, $region_end, $info) = (split /\t/)[0,1,2,3];
	
	my $cur_peak = {
    	chrom => $cur_chrom,
    	begin => $region_begin,
    	end => $region_end,
    	entry => $info				
	};
	$peaks[$peaks_counter] = $cur_peak;
	$peaks_counter++;
}
close calls_file;
printf "Read %d peaks\n", scalar(@peaks);

my $peaks_output_fname = $output_fname_prefix . ".peaks.". $window_size. "bp.fa";
if(-e $peaks_output_fname){
    unlink($peaks_output_fname);
}

open peaks_output_file, "> $peaks_output_fname" || die "$peaks_output_fname: $!\n";

print "\n";
print "Processing calls.\n\n";

my $chrom_counter = -1;
my $cur_chrom_name = "";
my $cur_chrom = "";
my $saved_peaks = 0;
my @peaks_mask;

for(my $i=0; $i<scalar(@peaks); $i++){
    $peaks_mask[$i] = 0;    
}

#maybe you could use bedtools to the following:
open ref_file, "< $ref_fname" || die "$ref_fname: $!\n";
while(<ref_file>){
    my $cur_line = $_;
    chomp($cur_line);
    
    if(length($cur_line) > 0){
		if(substr($cur_line, 0, 1) ne "#"){
	    	if(substr($cur_line, 0, 1) eq ">"){
				if($chrom_counter == -1){
				}
				else{
		    		# record the previous chromosome
		    		# process against the current chromosome
		    		printf "%.3f M bps\n", length($cur_chrom) / 1000000;
		    
		    		#$chrom_names[$chrom_counter] = $cur_chrom_name;
		    		#$chroms[$chrom_counter] = $cur_chrom;
		    		#print "\n";
				    for(my $i=0; $i<scalar(@peaks); $i++){
					my $this_peak_chrom = $peaks[$i]->{chrom};
					my $this_peak_begin = $peaks[$i]->{begin};
					my $this_peak_end = $peaks[$i]->{end};
					my $this_peak_entry = $peaks[$i]->{entry};
					
					if($this_peak_end - $this_peak_begin < $window_size){
					    my $diff = $window_size - ($this_peak_end - $this_peak_begin);
					    if($diff % 2 == 0){
							$this_peak_begin = $this_peak_begin - $diff/2;
							$this_peak_end = $this_peak_end + $diff/2;				
					    }
					    else{
							$this_peak_begin = $this_peak_begin - ($diff-1)/2;
							$this_peak_end = $this_peak_end + ($diff+1)/2;		
					    }
					}
					#my $cur_region_seq = substr($cur_chrom, $this_region_begin, $this_region_end);
					if($cur_chrom_name eq $this_peak_chrom){
					    	my $cur_peak_seq = substr($cur_chrom, $this_peak_begin, $this_peak_end - $this_peak_begin + 1);
					    	print peaks_output_file ">$this_peak_entry\n";
					    	print peaks_output_file "$cur_peak_seq\n";
					    	$saved_peaks++;
					    	#print "\rpeaks saved: $saved_peaks / $i";
					    	$|++;
					    	$peaks_mask[$i] = 1;
						}
		    		}
				}
				$cur_chrom_name = substr($cur_line, 1, length($cur_line) - 1);
				$cur_chrom = "";
				$chrom_counter++;	
				
				print "reading $cur_chrom_name : ";
				$|++;
	    	}
		    else{
				if($chrom_counter >= 0){
			    	$cur_chrom = $cur_chrom . $cur_line;
				}
		    }
		}
		else{
	    	print "skipping comment\n";
		}
    }
    else{
		print "skipping empty line\n";
    }
}

printf "%.3f M bps\n", length($cur_chrom) / 1000000;
print "peaks_saved: $saved_peaks\n";


print "\n";
for(my $i=0; $i<scalar(@peaks); $i++){
    my $this_peak_chrom = $peaks[$i]->{chrom};
    my $this_peak_begin = $peaks[$i]->{begin};
    my $this_peak_end   = $peaks[$i]->{end};
    my $this_peak_entry = $peaks[$i]->{entry};
    if($this_peak_end - $this_peak_begin < $window_size){
		my $diff = $window_size - ($this_peak_end - $this_peak_begin);
		if($diff % 2 == 0){
	    	$this_peak_begin = $this_peak_begin - $diff/2;
	    	$this_peak_end = $this_peak_end + $diff/2;				
		}
		else{
	  		$this_peak_begin = $this_peak_begin - ($diff-1)/2;
	   		$this_peak_end = $this_peak_end + ($diff+1)/2;		
		}
	}
    
    #my $cur_region_seq = substr($cur_chrom, $this_region_begin, $this_region_end);
    if($cur_chrom_name eq $this_peak_chrom){
	my $cur_peak_seq = substr($cur_chrom, $this_peak_begin, $this_peak_end - $this_peak_begin + 1);
	print peaks_output_file ">$this_peak_entry\n";
	print peaks_output_file "$cur_peak_seq\n";
	$saved_peaks++;
	print "\rpeaks saved: $saved_peaks / $i";
	$|++;
	$peaks_mask[$i] = 1;
    }
}
close ref_file;

close peaks_output_file;
#--------------------------------------------



#
#
##this concept of region has to disappear
#my @regions;
#my $regions_counter = 0;
#
#my @peaks;
#my $peaks_counter = 0;
#
#while( <calls_file> ){
#    chomp;
#    my $cur_line = $_;
#    if(length($cur_line) > 0){
#	if(substr($cur_line, 0, 1) eq "R"){
#	    my @cur_region_fields = split(/ /, $cur_line);
#	    if(scalar(@cur_region_fields) >= 3){
#		my $cur_chrom = $cur_region_fields[1];
#		my $cur_range = $cur_region_fields[2];
#		my $region_begin = -1;
#		my $region_end = -1;
#
#		my $updated_entry = "";
#		for(my $i=0; $i<length($cur_line); $i++){
#		    if(substr($cur_line, $i, 1) eq " "){
#			$updated_entry = $updated_entry . "~";
#		    }
#		    else{
#			$updated_entry = $updated_entry . substr($cur_line, $i, 1);
#		    }
#		}
#		
#		my @cur_range_fields = split(/-/, $cur_range);
#		if(scalar(@cur_range_fields) == 2){
#		    $region_begin = $cur_range_fields[0];
#		    $region_end = $cur_range_fields[1];
#
#		    if($region_end >= $region_begin && $region_begin >= 0){
#			my $cur_region = {
#			    chrom => $cur_chrom,
#			    begin => $region_begin,
#			    end => $region_end,
#			    entry => $updated_entry				
#			};
#			$regions[$regions_counter] = $cur_region;
#			$regions_counter++;
#		    }
#		    else{
#			print "error in the region range\n";
#
#			exit(1);
#		    }
#		}
#		else{
#		    print "Error in the region range.\n";
#		    print "Expected 2 fields but found $cur_range\n";
#		    exit(1);
#		}
#	    }
#	}
#	elsif(substr($cur_line, 0, 1)  eq "P"){
#	    my @cur_peak_fields = split(/ /, $cur_line);
#	    if(scalar(@cur_peak_fields) >= 3){
#		my $cur_chrom = $cur_peak_fields[1];
#		my $cur_coord = $cur_peak_fields[2];
#		
#		my $updated_entry = "";
#		for(my $i=0; $i<length($cur_line); $i++){
#		    if(substr($cur_line, $i, 1) eq " "){
#			$updated_entry = $updated_entry . "~";
#		    }
#		    else{
#			$updated_entry = $updated_entry . substr($cur_line, $i, 1);
#		    }
#		}
#		
#		my $cur_peak = {
#		    chrom => $cur_chrom,
#		    coord => $cur_coord,
#		    entry => $updated_entry
#		};
#
#		$peaks[$peaks_counter] = $cur_peak;
#		$peaks_counter++;
#	    }
#	}
#    }
#}
#
#close calls_file;
#printf "Read %d regions\n", scalar(@regions);
#printf "Read %d peaks\n", scalar(@peaks);
#
#my $peaks_output_fname = $output_fname_prefix . ".peaks.". $window_size. "bp.fa";
#my $regions_output_fname = $output_fname_prefix . ".regions." . $window_size . "bp.fa";
#
#if(-e $peaks_output_fname){
#    unlink($peaks_output_fname);
#}
#if(-e $regions_output_fname){
#    unlink($regions_output_fname);
#}
#
#open peaks_output_file, "> $peaks_output_fname" || die "$peaks_output_fname: $!\n";
#open regions_output_file, "> $regions_output_fname" || die "$regions_output_fname: $!\n";
#
##print "Reading the reference genome.\n\n";
#print "\n";
#print "Processing calls.\n\n";
#
#
##my @chroms;
##my @chrom_names;
#open ref_file, "< $ref_fname" || die "$ref_fname: $!\n";
#
#my $chrom_counter = -1;
#my $cur_chrom_name = "";
#my $cur_chrom = "";
#
#my $saved_regions = 0;
#my $saved_peaks = 0;
#
#my @peaks_mask;
#my @regions_mask;
#
#for(my $i=0; $i<scalar(@peaks); $i++){
#    $peaks_mask[$i] = 0;    
#}
#for(my $i=0; $i<scalar(@regions); $i++){
#    $regions_mask[$i] = 0;
#}
#
#
#
#while(<ref_file>){
#    my $cur_line = $_;
#    chomp($cur_line);
#    
#    if(length($cur_line) > 0){
#	if(substr($cur_line, 0, 1) ne "#"){
#	    if(substr($cur_line, 0, 1) eq ">"){
#		
#		
#		if($chrom_counter == -1){
#		    
#		}
#		else{
#		    # record the previous chromosome
#		    # process against the current chromosome
#		    
#		    printf "%.3f M bps\n", length($cur_chrom) / 1000000;
#		    
#		    #$chrom_names[$chrom_counter] = $cur_chrom_name;
#		    #$chroms[$chrom_counter] = $cur_chrom;
#		    #print "\n";
#		    for(my $i=0; $i<scalar(@regions); $i++){
#			my $this_region_chrom = $regions[$i]->{chrom};
#			my $this_region_begin = $regions[$i]->{begin};
#			my $this_region_end = $regions[$i]->{end};
#			my $this_region_entry = $regions[$i]->{entry};
#			if($this_region_end - $this_region_begin < $window_size){
#			    my $diff = $window_size - ($this_region_end - $this_region_begin);
#			    if($diff % 2 == 0){
#				$this_region_begin = $this_region_begin - $diff/2;
#				$this_region_end = $this_region_end + $diff/2;				
#			    }
#			    else{
#				$this_region_begin = $this_region_begin - ($diff-1)/2;
#				$this_region_end = $this_region_end + ($diff+1)/2;		
#			    }
#			}
#
#			#my $cur_region_seq = substr($cur_chrom, $this_region_begin, $this_region_end);
#			if($cur_chrom_name eq $this_region_chrom){
#			    my $cur_region_seq = substr($cur_chrom, $this_region_begin, $this_region_end - $this_region_begin + 1);
#			    print regions_output_file ">$this_region_entry\n";
#			    print regions_output_file "$cur_region_seq\n";
#			    $saved_regions++;
#			    #print "\rregions saved: $saved_regions / $i";
#			    $|++;
#			    $regions_mask[$i] = 1;
#			}
#		    }
#
#		    #print "\n";
#
#		    for(my $i=0; $i<scalar(@peaks); $i++){
#			my $this_peak_chrom = $peaks[$i]->{chrom};
#			my $this_peak_coord = $peaks[$i]->{coord};
#			my $this_peak_entry = $peaks[$i]->{entry};
#
#			my $region_begin = $this_peak_coord - $window_size / 2;
#			my $region_end = $this_peak_coord + $window_size / 2;
#
#			if($cur_chrom_name eq $this_peak_chrom){
#			    my $cur_region_seq = substr($cur_chrom, $region_begin, $region_end - $region_begin + 1);
#			    print peaks_output_file ">$this_peak_entry\n";
#			    print peaks_output_file "$cur_region_seq\n";
#			    $saved_peaks++;
#			    #print "\rpeaks saved: $saved_peaks / $i";
#			    $|++;
#			    $peaks_mask[$i] = 1;
#			}
#		    }
#		    #print "\n";
#		    #print "saved $saved_regions regions\n";
#		    #print "\n\n";
#		}
#
#
#		$cur_chrom_name = substr($cur_line, 1, length($cur_line) - 1);
#		$cur_chrom = "";
#		$chrom_counter++;	
#		
#		print "reading $cur_chrom_name : ";
#		$|++;
#	    }
#	    else{
#		if($chrom_counter >= 0){
#		    $cur_chrom = $cur_chrom . $cur_line;
#		}
#	    }
#	}
#	else{
#	    print "skipping comment\n";
#	}
#    }
#    else{
#	print "skipping empty line\n";
#    }
#}
#
#
#
#
#
#
#
#printf "%.3f M bps\n", length($cur_chrom) / 1000000;
#
#print "regions_saved: $saved_regions\n";
#print "peaks_saved: $saved_peaks\n";
#
#
#
##$chrom_names[$chrom_counter] = $cur_chrom_name;
##$chroms[$chrom_counter] = $cur_chrom;
#print "\n";
#for(my $i=0; $i<scalar(@regions); $i++){
#    my $this_region_chrom = $regions[$i]->{chrom};
#    my $this_region_begin = $regions[$i]->{begin};
#    my $this_region_end = $regions[$i]->{end};
#    my $this_region_entry = $regions[$i]->{entry};
#    if($this_region_end - $this_region_begin < $window_size){
#	my $diff = $window_size - ($this_region_end - $this_region_begin);
#	if($diff % 2 == 0){
#	    $this_region_begin = $this_region_begin - $diff/2;
#	    $this_region_end = $this_region_end + $diff/2;				
#	}
#	else{
#	    $this_region_begin = $this_region_begin - ($diff-1)/2;
#	    $this_region_end = $this_region_end + ($diff+1)/2;		
#	}
#    }
#    
#    #my $cur_region_seq = substr($cur_chrom, $this_region_begin, $this_region_end);
#    if($cur_chrom_name eq $this_region_chrom){
#	my $cur_region_seq = substr($cur_chrom, $this_region_begin, $this_region_end - $this_region_begin + 1);
#	print regions_output_file ">$this_region_entry\n";
#	print regions_output_file "$cur_region_seq\n";
#	$saved_regions++;
#	print "\rregions saved: $saved_regions / $i";
#	$|++;
#	$regions_mask[$i] = 1;
#    }
#}
#
##-------------------------------------------
#
#print "\n";
#
#for(my $i=0; $i<scalar(@peaks); $i++){
#    my $this_peak_chrom = $peaks[$i]->{chrom};
#    my $this_peak_coord = $peaks[$i]->{coord};
#    my $this_peak_entry = $peaks[$i]->{entry};
#    
#    my $region_begin = $this_peak_coord - $window_size / 2;
#    my $region_end = $this_peak_coord + $window_size / 2;
#    
#    if($cur_chrom_name eq $this_peak_chrom){
#	my $cur_region_seq = substr($cur_chrom, $region_begin, $region_end - $region_begin + 1);
#	print peaks_output_file ">$this_peak_entry\n";
#	print peaks_output_file "$cur_region_seq\n";
#	$saved_peaks++;
#	print "\rpeaks saved: $saved_peaks / $i";
#	$|++;
#	$peaks_mask[$i] = 1;
#    }
#}
##print "\n";
##print "saved $saved_regions regions\n";
#print "\n\n";
#
#close ref_file;
#
#close peaks_output_file;
#close regions_output_file;
##-------------------------------------------
#print "\n";
#
#my $regions_skipped = "false";
#my $peaks_skipped = "false";
#
#for( my $i=0; $i<scalar(@regions); $i++){
#    if($regions_mask[$i] != 1){
#	$regions_skipped = "true";
#    }
#}
#
#for(my $i=0; $i<scalar(@peaks); $i++){
#    if($peaks_mask[$i] != 1){
#	$peaks_skipped = "true";
#    }
#}
#
#if($regions_skipped eq "true"){
#    print "There were skipped regions.\n";
#}
#else{
#    print "All regions were saved.\n";
#}
#
#if($peaks_skipped eq "true"){
#    print "There were skipped peaks.\n";
#}
#else{
#    print "All peaks were saved.\n";
#}
#
##printf "\n";
##printf "read %d chromosomes\n", scalar(@chroms);
##print "\n";
#
##if(-e $accept_fname){ unlink("$accept_fname");}
##open accept_file, "> $accept_fname" || die "$accept_fname: $!\n";
#
##open regions_file, "< $regions_fname" || die "$regions_fname: $!\n";
##}
#
#
##close regions_file;
##close accept_file;
##45close reject_file;
