#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;
use List::Util qw(min max sum reduce);

#Get the .rData file you got from do_funseq_collect_GWAS_effect_dir.pl and postprocess according to the following filters

#1-remove non-european or yri "INIT_SAMPLE" AND "REP_SAMPLE"
#2-remove cases with "or" not reported 
#3-remove pvals < 0.00000005
#4-remove terms which are not enriched?
#5-remove variants in the MHC
#6-leave only 1 variant per LD block?

my $PVAL_THRS = '0.00000005';

#mhc exclusion-------------
my $MHC_chr = '6';
my $MHC_start = '29540169';
my $MHC_end = '33215544';
#mhc exclusion-------------

my $INPUT_GWAScat_LDmapped_LD1 = "/net/isi-scratch/giuseppe/indexes/GWAS_Gwascat/UCSC_gwas_catalog/gwasCatalog_UCSC_mapped_LDblocks_BEAGLE_CEU_r2_1.0_b37.tsv.gz";
my $INPUT_GWAScat_LDmapped_LD08 = "/net/isi-scratch/giuseppe/indexes/GWAS_Gwascat/UCSC_gwas_catalog/gwasCatalog_UCSC_mapped_LDblocks_BEAGLE_CEU_r2_0.8_b37.tsv.gz";

my $INPUT_GWAScat_LDmapped = $INPUT_GWAScat_LDmapped_LD1;

my $INPUT_ALLELESEQ_RAW = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/_RAW/interestingHets_vdrchipexo.txt"; #b37
my $BEDTOOLS = "/net/isi-scratch/giuseppe/tools/bedtools-2.22.1/bin/bedtools";

my $infile_rdata;
my $infile_terms;
my $LDCLUMP;

GetOptions(
        'i=s'     =>\$infile_rdata,
        't=s'     =>\$infile_terms,
        'ldclump' =>\$LDCLUMP
);

#$infile_rdata = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/FIGURE6/FIG6_Output_noDBRECUR_ALL.vcf_gwas_effectLD1.Rdata";
#$infile_terms = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/FIGURE6/enriched_terms_all_only_MS_C_IBD_T1D.txt";

if(!$infile_rdata){
     print "USAGE: do_funseq_collect_GWAS_effect_dir_postproc.pl -i=<IN_RDATA> -t=<IN_TERMS> -ldclump\n";
     print "<IN_RDATA> tsv file, output of do_funseq_collect_GWAS_effect_dir.pl\n";
     print "(optional)<IN_TERMS> list of terms to keep because the VDR-BVs are ENRICHED in them\n";
     print "(optional)<ldclump> if set, only one VDR-BV per LD block (r2=1.0,BROAD) will be kept\n";
     exit 0;
}

my($basename, $directory) = fileparse($infile_rdata);
$basename =~ s/(.*)\..*/$1/;

#temp bed file to put filtered variants in a bed and select by LD
my $out_temp_bed        = $directory . 'TEMP_'      . $basename . '_list.bed';
my $out_temp_buffer     = $directory . 'TEMP_'       . $basename . '_buffer.bed';
my $output_full;
my $output_ldclumped;
if($infile_terms){
	$output_full         = $directory .  'TABLE_pp_' . $basename . '_nonLDclumped_onlyenrichterms.rData'; 
	$output_ldclumped    = $directory .  'TABLE_pp_' . $basename . '_LDclumped_onlyenrichterms.rData'; 		
}else{
	$output_full         = $directory .  'TABLE_pp_' . $basename . '_nonLDclumped.rData'; 
	$output_ldclumped    = $directory .  'TABLE_pp_' . $basename . '_LDclumped.rData'; 	
}

#GLOBALS
my %seen;
my %data;
my %genomic_positions; #get bed from this
my %data_nold;
my %data_trait2ld2var;
my %outdump;
my %candidate_vars;

my @ENRICHED_TRAITS;
if($infile_terms){
	open (my $instream,  q{<}, $infile_terms) or die("Unable to open $infile_terms : $!");
	while(<$instream>){
		chomp;
		push(@ENRICHED_TRAITS, $_);
	}
	close $instream;	
}

#input fields
#0 SNP
#1 Chr
#2 pos
#3 MainGroup:pval
#4 MainGroup:or
#5 MainGroup:lower_ci
#6 MainGroup:upper_ci
#7 MainGroup:caf
#8 MainGroup:CEU_DAF
#9 CI_DATA
#10 BINDING_MAGNITUDE
#11 SAMPLE_ID
#12 ASB_PVAL
#13 BINDING_DIR
#14 RISK_DIR
#15 INIT_SAMPLE
#16 REP_SAMPLE
#17 EFFECT
#18 TRAIT_DISEASE
#19 AUTHOR
#20 RC
my $HEADER;
open (my $instream,  q{<}, $infile_rdata) or die("Unable to open $infile_rdata : $!");
while(<$instream>){
	if($_ =~ /^SNP/){
		$HEADER = $_;
		next;
	}
	chomp;
	my @fields = split("\t",$_);
	my $coords = $fields[1] . '-' . $fields[2];
	
	#check fields
	#test 0 MHC-----------------------
	if($fields[1] eq $MHC_chr){
		next if( ($fields[2] > $MHC_start) &&  ($fields[2] < $MHC_end) );
	}
	#test 1: pval---------------------
	next if(!$fields[3]);
	next unless($fields[3] <= $PVAL_THRS);
	#test 2: OR-----------------------
	next if(!$fields[4] || $fields[4] eq ''  || $fields[4] =~ /NR/);
	#test SAMPLE non reported---------
	next if($fields[15] =~ /NA/ || $fields[16] =~ /NA/);
	#test INIT SAMPLE chinese (add other cases?)
	#next if($fields[13] =~ /chinese/i || $fields[13] =~ /korea/i || $fields[13] =~ /japan/i || $fields[13] =~ /mexican/i || $fields[13] =~ /asian/i);
	#next if($fields[14] =~ /chinese/i || $fields[14] =~ /korea/i || $fields[14] =~ /japan/i || $fields[14] =~ /mexican/i || $fields[14] =~ /asian/i);
	next unless($fields[15] =~ /european/i);
	#test trait is amongst enriched list-------------------
	if($infile_terms){
		my $found;
		foreach my $enriched_term (@ENRICHED_TRAITS){
			if($fields[18] =~ qr/$enriched_term/i){
				$found = 1;
			}
		}
		next unless $found;
	}
	
	#save coords
	$genomic_positions{$coords} = 1;
	#save rsID with counts
	$seen{$fields[0]}++;
	#save full data
	$data{$_} = 1;
}
close $instream;

##########################
#0 save non LD-clumped out
##########################
foreach my $lines (keys %data){
	my @fields = split("\t",$lines);
	my $rsid = shift @fields;
	
	if($seen{$rsid} > 1){
		$rsid = $rsid . '_' . $seen{$rsid}--;
	}else{
		$rsid = $rsid . '_1';
	}
	my $newline = $rsid . "\t" . join("\t",@fields);
	#print $outstream $rsid . "\t" . join("\t",@fields) . "\n";
	$outdump{$newline} = 1;
}

open (my $outstream,  q{>}, $output_full) or die("Unable to open $output_full: $!");
print $outstream $HEADER;
foreach my $item (keys %outdump){
	print $outstream $item, "\n";
}
%outdump = ();
close $outstream;

my $tot_cand = scalar (keys %seen);
print STDOUT "Before LD-clumping, $tot_cand rsIDs are retained.\n";

exit unless $LDCLUMP;





#LD clumping
##########################
#1 create bed of positions
##########################
open ($outstream,  q{>}, $out_temp_bed) or die("Unable to open $out_temp_bed: $!");		
foreach my $item (sort keys %genomic_positions){
	my ($chr, $pos) = split("-", $item);
	print $outstream   $chr . "\t" . ($pos - 1) . "\t" . $pos . "\n"; 
}
close $outstream;


################
#2 create a hash to get all the fg/bg variants which are NOT in an LD block
################
system "$BEDTOOLS intersect -c -a $out_temp_bed -b $INPUT_GWAScat_LDmapped > $out_temp_buffer";
store_variants_not_in_LDblock($out_temp_buffer, \%data_nold);
unlink $out_temp_buffer;
###################
#3 for each trait allowed, map LD block with at least a fg/bg variant in them to variant coordinate
###################
system "$BEDTOOLS intersect -wo -a $INPUT_GWAScat_LDmapped -b $out_temp_bed > $out_temp_buffer";
map_LDblocks_to_variants($out_temp_buffer, \%data_trait2ld2var);

###################
#4 Get per-sample affinity binding phenotype info (to find the best VDR-BV in an LD block)
###################
my %position2sample2readdepth;
get_binding_affinity_levels(\%position2sample2readdepth);

####################
#5 go through all the LD blocks intersecting vdr-bvs and choose 1 foreground var for each
####################
foreach my $trait (sort keys %data_trait2ld2var){
	my $found;
	#check trait is allowed, otherwise skip
	if($infile_terms){
		foreach my $enriched_term (@ENRICHED_TRAITS){
			if($trait =~ qr/$enriched_term/i){
				$found = 1;
			}
		}
	}
	if(!$found){
		print STDOUT "LDclumping: skipping trait: $trait because it's not in the enriched list..\n";
		next;
	}
	
	foreach my $ld_block (keys %{ $data_trait2ld2var{$trait} }){
		my @vars_in_ld; 
		my $candidate_var;	
		
		foreach my $var_in_ld (sort keys %{ $data_trait2ld2var{$trait}{$ld_block} }  ){
			push(@vars_in_ld, $var_in_ld);
		}	
		if(@vars_in_ld == 1){
			$candidate_var = $vars_in_ld[0];
		}else{
			#call routine to select the "best" and fill candidate_vdrucbv
			$candidate_var = get_best_var(@vars_in_ld);
		}	
		
		if($candidate_var){
			$candidate_vars{$candidate_var} = 1;
		}else{
			print STDOUT "WARNING: candidate variant undefined. Skipping..\n";
			next;
		}		
		
	}
}

###############
#6 save the foreground and background variables which are not in LD blocks
###############
foreach my $item (sort keys %data_nold){
	$candidate_vars{$item} = 1;
}
#############
#7 print foreground and background
#############
$tot_cand = scalar keys %candidate_vars;
print STDOUT "After LD-clumping, $tot_cand variant positions are retained.\n";

##############
#8 relabel rsIDs which appear more than 1 times (the graph program at http://visualization.ritchielab.psu.edu/phenograms/plot will not plot them otherwise!!)
##############
%seen = ();
foreach my $line (keys %data){
	my @fields = split("\t",$line);
	my $coord = $fields[1] . '-' . $fields[2];
	$seen{$fields[0]}++ if($candidate_vars{$coord});
}
$tot_cand = scalar (keys %seen);
print STDOUT "After LD-clumping, $tot_cand rsIDs are retained.\n";
###############
#9 print output
###############
foreach my $line (keys %data){
	my @fields = split("\t",$line);
	my $coord = $fields[1] . '-' . $fields[2];

	if($candidate_vars{$coord}){
		my $rsid = shift @fields;
		
		if($seen{$rsid} > 1){
			$rsid = $rsid . '_' . $seen{$rsid}--;
		}else{
			$rsid = $rsid . '_1';
		}	
		my $line = $rsid . "\t" . join("\t",@fields);
		#print $outstream $rsid . "\t" . join("\t",@fields) . "\n";
		$outdump{$line} = 1;
		
	}
	#my $rsid = shift @fields;
}
open ($outstream,  q{>}, $output_ldclumped) or die("Unable to open $output_ldclumped: $!");
print $outstream $HEADER;
foreach my $item (keys %outdump){
	print $outstream $item, "\n";
}

close $outstream;
unlink $out_temp_buffer;
unlink $out_temp_bed;

#################
#subs
#################
sub store_variants_not_in_LDblock{
	my ($file, $hash) = @_;
	
	#format:
	#22	27068481	27068482	1
	#22	28807624	28807625	1
	#22	37258266	37258267	1
	#22	37258266	37258267	0  <- I want a hash with these
	open (my $instream,  q{<}, $file) or die("Unable to open $file: $!");		
	while(<$instream>){
		chomp;
		next if($_ eq '');	
		my ($chr, $pos, $in_ld) = (split /\t/)[0,2,3];
		if($in_ld == 0){
			my $coord = $chr . '-' . $pos;
			$$hash{$coord} = 1;
			next;
		}
	}
	close $instream;
	return;	
}

sub map_LDblocks_to_variants{
	my ($file, $hash) = @_;
	
	#format:
	#chr	ld_start	ld_end	info_snps	vdrhcbv_chr	vdrhcbv_start	vdrhcbv_stop	n_bases
	#1	
	#197611681
	#197631141
	#197631141__rs2488389__23128233__Jostins L__2012-11-01__Nature__Host-microbe interactions have shaped the genetic architecture of inflammatory bowel disease.__Inflammatory bowel disease__12,924 European ancestry cases, 21,442 European ancestry controls__25,683 European ancestry cases, 17,015 European ancestry controls__1q31.3__C1orf53__rs2488389-A__0.22__8E-13____1.12__[1.077-1.153]__Affymetrix & Illumina [1.23 million](imputed)__N
	#1
	#197619067
	#197619068
	#1
	open (my $instream,  q{<}, $file) or die("Unable to open $file: $!");		
	while(<$instream>){
		chomp;
		next if($_ eq '');	
		my ($ld_chr, $ld_start, $ld_stop, $gwas_data, $var_chr, $var_stop) = (split /\t/)[0,1,2,3,4,6];
		unless($ld_chr eq $var_chr){
			print STDERR "map_LDblocks_to_variants() - Error: chromosomes do not coincide: $ld_chr vs $var_chr. Aborting..\n";
			return undef;
		}
		
		#get the term
		my @gwas_fields = split('__', $gwas_data);
		my $this_trait = $gwas_fields[7];
		
#		if($infile_terms){
#			my $found;
#			foreach my $enriched_term (@ENRICHED_TRAITS){
#				if($this_trait =~ qr/$enriched_term/i){
#					$found = 1;
#				}
#			}
#			if(!$found){
#				print STDERR "It should not be here. This trait is $this_trait. It's not amongst the enriched traits. Aborting..\n";
#				exit -1;
#			}
#		}
		my $ld_coord = $ld_chr . '-' . $ld_start . '-' . $ld_stop;
		my $var_coord = $ld_chr . '-' . $var_stop;
		$$hash{$this_trait}{$ld_coord}{$var_coord} = 1;
	}
	close $instream;
	return;
}


#chr-position => samplename => alleles / read counts info
sub get_binding_affinity_levels{
	my ($hash) = @_;
	
	open (my $instream,  q{<}, $INPUT_ALLELESEQ_RAW) or die("Unable to open $INPUT_ALLELESEQ_RAW : $!");
	
	while(<$instream>){
		chomp;
		next if($_ =~ /^sample/); #header
		next if($_ eq '');
		
		#format
		#sample	chrm	snppos	ref	mat_gtyp	pat_gtyp	c_gtyp	phase	mat_all	pat_all	cA	cC	cG	cT	winning	SymCls	SymPval	BindingSite	cnv
		#NA06986	1	1080920	G	S	W	R	PHASED	G	A	2	0	7	0	M	Sym	0.1796875	1	1.0
		my ($sample_id, $chr,$pos, $ref, $m_allele, $p_allele, $cA, $cC, $cG, $cT, $winner, $symcls) = (split /\t/)[0,1,2,3,8,9,10,11,12,13,14,15];
		
		next if(!$chr);
		next if(!$pos);
		next unless ($symcls =~ /Asym/);
		#$chr = 'chr' . $chr;
		my $coordinate_id =  $chr . '-' . $pos; #b37
		
		my $read_data = join(",", $ref, $cA, $cC, $cG, $cT);
		$$hash{$coordinate_id}{$sample_id} = $read_data;
	}
	close $instream;
	return;
}


sub get_best_var{
	my (@ld_variant_vector) = @_;
	my %vdrbv_to_magnitude_hash;
	
	#need to cycle by sample
	foreach my $vdrbv (@ld_variant_vector){
		my @persample_dba_magns;
		#print "\n\nVDRBV: " . $vdrbv . "\n";
		
		#check if VDR-BV is VDR-QTL (it might not be in the ASB hash)
		unless(defined $position2sample2readdepth{$vdrbv}){
			#then this is a QTL. There should be only a dozen or so in the VDR-hcBV set. Keep it?
			print STDERR "WARNING: get_best_var(): VDR-BV: $vdrbv is not in the VDR-ASB set so it's a VDR-QTL.\n";
			print STDERR "Coverage unknown: assigning magnitude 0.\n";
			$vdrbv_to_magnitude_hash{$vdrbv} = 0;
			next;
		}
		
		foreach my $sample (sort keys %{ $position2sample2readdepth{$vdrbv} }  ){
				my ($ref, $cA,$cC,$cG,$cT) = split(",", $position2sample2readdepth{$vdrbv}{$sample});
				if( ($cA eq 0) && ($cC eq 0) && ($cG eq 0) && ($cT eq 0) ){
					print STDERR "get_vdrucbv(): error: 0 read coverage for all possible nucleotides at this event ($cA,$cC,$cG,$cT). Aborting..\n";
					exit -1;
				}
				my @allele_rd = ($cA, $cC, $cG, $cT);
				my @sorted_allele_rd = sort {$b <=> $a}  @allele_rd;
				#print 'sample: ' . $sample . ': ' . join(",",@sorted_allele_rd) . "\n";
				
				#pop the largest
				my $max_coverage_allele = shift(@sorted_allele_rd);
				if($max_coverage_allele == 0){
					print STDERR "WARNING: get_best_var(): max coverage allele is 0 in this sample. Skipping..\n";
					next;
				}
				#pop the second largest
				my $second_coverage_allele = shift(@sorted_allele_rd);
				
				#if a third is non zero, skip this sample altogether (mismapping)
#				my $MISMATCH_FLAG;
#				foreach my $item (@sorted_allele_rd){
#					if($item != 0){
#						print STDERR "WARNING: get_vdrucbv(): more than 2 alleles have coverage: $item. Skipping this sample..\n";
#						$MISMATCH_FLAG = 1;
#					}
#				}
#				next if($MISMATCH_FLAG);
				
				my $dba_magnitude_thissample = $max_coverage_allele - $second_coverage_allele;
				push (@persample_dba_magns, $dba_magnitude_thissample);
		}
		#print 'Vector of per-sample magnitudes: ' .  join(',', @persample_dba_magns) . "\n";
		if(@persample_dba_magns == 1){
			$vdrbv_to_magnitude_hash{$vdrbv} = $persample_dba_magns[0];
		}else{
			my $dba_magnitude = mean(@persample_dba_magns);
			#print 'Mean per-sample magnitude: ' .  $dba_magnitude, "\n";
			$vdrbv_to_magnitude_hash{$vdrbv} = $dba_magnitude;
		}
		
	}
	#return the hash key with the largest value;
	return reduce { $vdrbv_to_magnitude_hash{$a} > $vdrbv_to_magnitude_hash{$b} ? $a : $b } keys %vdrbv_to_magnitude_hash;
	#print '------------VDR-hucBV: ' . $out . "\n";
}

#------mean
sub mean {
    return sum(@_)/@_;
}
