#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;

#Chris wants me to cross reference VDR-BVs and eQTLs from the geuvadis data:
#I think in the abstract we can say that we predict XX genes to be differentially under the control of VDR in the human population (i.e. VDR-BVs coincident with LCL #eQTLs, despite being calcitriol-minus) and YY genes that are both under the control of VDR and contribute to disease susceptibility (i.e. VDR-BVs coincident with both #LCL eQTLs and GWAS SNPs [or in LD]). We're going to have a hard time over this with Ram, but the data justify this, IMHO.
#Re: LCL eQTLs
#Are there not tables of LCL eQTL SNPs cross-referenced with the genes whose expression they correlate with? If so, this is what I am interested in. I don't see why you #should recall eQTLs using LCL RNA-Seq data - surely this has already been done!?

#The sign denotes the direction of the nonreference allele (i.e. rvalue<0 means that nonreference allele has lower expression) 
#CONVERT THIS TO ANCESTRAL/DERIVED

#inputs
#1 - output.tsv from funseq2, indicates the positions (NO DBRECUR if you want only VDR-hcBVs)
#2 - 1kg data, gives the ancestral -> derived
#3 - output from alleleseq, gives the binding direction
#4 - GEUVADIS file

#what to do when x samples are 'lob' and one is 'gob'? ATM I discard them, unless the x >= 4y (take x) or y >= 4x (take y)

#METHOD
#VDR-QTL-GENE relationships
#1:m - VDR-BV is EQTL for m genes - test and report concordance for all of them
#m:1 - Gene is affected by multiple VDR-BV/QTL: choose the most significant
#the "best" datasets should already contain the most significant variant for each gene
#"best association per each gene (for exon, transcript ratio) or unit (gene, repeat, miRNA). If there are several best associating variants with the same p-value, one of them has been chosen randomly. "

#to calculate a fdr, get 1000 random sets of NON-VDR-BVs, same number of VDR-BVs intersecting the eQTL, and compute a distribution of the overlap
my $INPUT_VARIANTS = '/net/isi-scratch/giuseppe/VDR/VARIANTS/ALLELESEQ_20SAMPLES/ALLELESEQ_20samples_merged.vcf.gz'; #b37
my $INPUT_ALLELESEQ_RAW = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/RESULTS/_RAW/interestingHets_vdrchipexo.txt"; #b37


#typical input structure
#SNP_ID
#ID
#GENE_ID
#PROBE_ID
#CHR_SNP
#CHR_GENE
#SNPpos
#TSSpos
#distance
#rvalue
#pvalue
#log10pvalue

#rs3832000
#-
#ENSG00000142632.10
#ENSG00000142632.10
#1
#1
#16533832
#16539140
#5308
#-0.765767800472289
#3.68318824172032e-18
#17.433776084498

#1	SNP_ID : Variant identifier according to dbSNP137; position-based identifier for variants that are not in dbSNP (see Supplementary material pp 45) 
#2	ID : Null (-)
#3	GENE_ID : Gene identifier according to Gencode v12, miRBase v18, repeats based on their start site
#4	PROBE_ID : Quantitative trait identifier; the same as GENE_ID expect for:
#		Exons: GENEID_ExonStartPosition_ExonEndPosition
#		Transcript ratios: Transcript identifier according to Gencode v12
#5	CHR_SNP : Chromosome of the variant
#6	CHR_GENE : Chromosome of the quantitative trait
#7	SNPpos : Position of the variant
#8	TSSpos : Transcription start site of the gene/QT
#9	Distance : | SNPpos - TSSpos | 
#10	rvalue	: Spearman rank correlation rho (calculated from linear regression slope). The sign denotes the direction of the nonreference allele (i.e. rvalue<0 means that nonreference allele has lower expression) 
#11	pvalue : Association p-value
#12	log10pvalue : -log10 of pvalue

#loss of binding, gain of binding
my $LOB = 'lob'; #ancestral -> alternative reduces binding
my $GOB = 'gob'; #ancestral -> alternative increase binding
#loss of expression, gain of expression
my $LOE = 'loe'; #ancestral -> alternative reduces expression  
my $GOE = 'goe'; #ancestral -> alternative increases expression
my $ACTIVATE = 'VDR_ACTIVATES';
my $REPRESS = 'VDR_REPRESSES';
#LOB + LOE -> transition from ancestral to alternative cause:  1) decrease of VDR binding 2) decrease of expression - VDR activator
#GOB & GOE -> transition from ancestral to alternative cause:  1) increase of VDR binding 2) increase of expression - VDR activator
#LOB & GOE -> transition from ancestral to alternative cause:  1) decrease of VDR binding 2) increase of expression - VDR repressor
#GOB & LOE -> transition from ancestral to alternative cause:  1) increase of VDR binding 2) decrease of expression - VDR repressor

my $input_funseq;
my $input_geuvadis;
my $sym_variants;			
my $motif_name = 'VDR_JASPAR';
GetOptions(
        'i_fun=s'  =>\$input_funseq,
        'i_geu=s'  =>\$input_geuvadis,
        'sym'      =>\$sym_variants
);
#allow to select items in class I sites only
#$input_funseq = "/net/isi-scratch/giuseppe/VDR/ALLELESEQ/funseq2/out_allsamples_plus_qtl_ancestral/Output_noDBRECUR.vcf";

#in figure s21b you show the intersection with ALL. However you have computed the intersection with the best as well.
#the following file cat'ed the best for gene, exon,and ttratio QTL
#$input_geuvadis = "/net/isi-scratch/giuseppe/indexes/GEUVADIS/BEST/d_raw/EUR373.alltypes.cis.FDR5.best.rs137.txt.gz";
#$input_geuvadis = "/net/isi-scratch/giuseppe/indexes/GEUVADIS/BEST/d_raw/YRI89.alltypes.cis.FDR5.best.rs137.txt.gz";
#all qtls
#$input_geuvadis = "/net/isi-scratch/giuseppe/indexes/GEUVADIS/d_raw/EUR373.alltypes.cis.FDR5.all.rs137.txt.gz";
#$input_geuvadis = "/net/isi-scratch/giuseppe/indexes/GEUVADIS/d_raw/YRI89.alltypes.cis.FDR5.all.rs137.txt.gz";

#$input_geuvadis = "/net/isi-scratch/giuseppe/indexes/GEUVADIS/d_raw/EUR373.gene.cis.FDR5.best.rs137.txt.gz";
#$input_geuvadis = "/net/isi-scratch/giuseppe/indexes/GEUVADIS/d_raw/YRI89.gene.cis.FDR5.best.rs137.txt.gz";
#$input_geuvadis = "/net/isi-scratch/giuseppe/indexes/GEUVADIS/d_raw/EUR373.gene.cis.FDR5.all.rs137.txt.gz";
#$input_geuvadis = "/net/isi-scratch/giuseppe/indexes/GEUVADIS/d_raw/YRI89.gene.cis.FDR5.all.rs137.txt.gz";

if(!$input_funseq){ #hg19(convert)
     print "USAGE: do_funseq_GEUVADIS_eQTL_direction.pl -i_fun=<INFILE_FUNSEQ> -i_geu=<INFILE_GEUVADIS> -sym\n";
     print "<INFILE_FUNSEQ> vcf output of funseq - NODBRECUR\n";
     print "<INFILE_GEUVADIS_BEST> raw tsv with expression data, best per gene\n";
     print "GEUVADIS files are in /net/isi-scratch/giuseppe/indexes/GEUVADIS\n";     
     print "(optional)<sym> carry out analysis for Output.vcf of non VDR-BVs variants\n";
     exit 1;
}
if(!$input_geuvadis){ #b37
     print "USAGE: do_funseq_GEUVADIS_eQTL_direction.pl -i_fun=<INFILE_FUNSEQ> -i_geu=<INFILE_GEUVADIS> -sym\n";
     print "<INFILE_FUNSEQ> vcf output of funseq - NODBRECUR\n";
     print "<INFILE_GEUVADIS_BEST> raw tsv with expression data, best per gene\n";
     print "GEUVADIS files are in /net/isi-scratch/giuseppe/indexes/GEUVADIS\n";
     print "(optional)<sym> carry out analysis for Output.vcf of non-VDR-BVs variants\n";
     exit 1;
}
#$infile = "/net/isi-scratch/giuseppe/tools/funseq2-1.0/out_allsamples_plus_qtl_ancestral/Output_noDBRECUR.vcf"; #hg19
#you only want the recurrent from here
my($basename_funseq, $directory_funseq) = fileparse($input_funseq);
my($basename_geuvadis, $directory_geuvadis) = fileparse($input_geuvadis);
$basename_funseq =~ s/(.*)\..*/$1/;
$basename_geuvadis =~ s/(.*)\..*/$1/;
my $output  = $directory_funseq . 'eQTL_GEUVADIS_' . $basename_funseq . '_' .  $basename_geuvadis . '.data';

###################
#1 get QTL data 
###################
#SNP_ID	ID	GENE_ID	PROBE_ID	CHR_SNP	CHR_GENE	SNPpos	TSSpos	distance	rvalue	pvalue	log10pvalue
#rs479844	-	ENSG00000254470.2	ENSG00000254470.2	11	11	65551957	65548273	3684	-0.238920740106044	3.07378298927364e-06	5.51232679717714
#rs479844	-	ENSG00000172818.5	ENSG00000172818.5	11	11	65551957	65554493	2536	-0.548355722485014	1.29308592493011e-30	29.8883726155039
my %snp2QTL_info;
tie *FILE,   'IO::Zlib', $input_geuvadis, "rb";
while (<FILE>)	{ 
	chomp;
	next if($_ =~ /^SNP_ID/);
	next if($_ eq '');

	my ($geu_snp_id, $geu_gene_id, $geu_snp_chr, $geu_snp_pos, $geu_rval, $geu_pval) = (split /\t/)[0,2,4,6,9,10];
	next if($geu_snp_id =~ /\./); #skip indels
	unless($geu_snp_id && $geu_gene_id && $geu_snp_chr && $geu_snp_pos && $geu_rval && $geu_pval){
		print STDERR "Warning, something is missing in this GEUVADIS line:\n";
		print STDERR $_, "\n";
		print STDERR "Skipping..\n";
		next;
	}
	my $geu_coordinate_id = $geu_snp_chr. '-' . $geu_snp_pos;
	#data structure
	$snp2QTL_info{$geu_coordinate_id}{$geu_snp_id}{$geu_gene_id}{RVAL} = $geu_rval;
	$snp2QTL_info{$geu_coordinate_id}{$geu_snp_id}{$geu_gene_id}{PVAL} = $geu_pval; #useless?
}
close FILE;


#####################
#2 -build hash to map chr-pos to ref, ancestral, alternate, and DAFs
#####################
#before saving the frequencies, you need to to KNOW if the ancestral allele is the ref or the alt
#if the ancestral allele is the ref, save the frequency as is
#if the ancestral allele is the alt, the frequency you have is for the ancestral. The derived will be (1 - freq)
my %variants_1kg;
tie *FILE,   'IO::Zlib', $INPUT_VARIANTS, "rb";
while (<FILE>)	{ 
	chomp;
	next if($_ eq '');
	next if($_ =~ /^#/);
	my $FLAG; # set to one if the ancestral is the alternate	
	
	my @fields = split("\t", $_);
	#skip indels
	next if(length($fields[3]) > 1);
	next if(length($fields[4]) > 1);
	
	my $key = $fields[0] . '-' . $fields[1];
	my $ref = uc($fields[3]); 
	my $alt = uc($fields[4]);
	my $info = $fields[7];
	
	#get ancestral allele info=====================
	my @info = split(";", $fields[7]);
	if($info[0] =~ /^AA=(.*)/){
		my $anc = uc($1);
		if($anc eq $alt){
			$variants_1kg{$key}{ANC} = $alt;
			$variants_1kg{$key}{DER} = $ref;
			$variants_1kg{$key}{SWAP} = 1;		
		}elsif($anc eq $ref){
			$variants_1kg{$key}{ANC} = $ref;
			$variants_1kg{$key}{DER} = $alt;			
		}else{
			next;
		}			
	}else{
		next;	
	}
}
close FILE;

###################
#4 get affinity binding phenotype info 
###################
#chr-position => samplename => alleles / read counts info
my %position2sample2readdepth;
open (my $instream,  q{<}, $INPUT_ALLELESEQ_RAW) or die("Unable to open $INPUT_ALLELESEQ_RAW : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^sample/); #header
	next if($_ eq '');

	#format
	#sample	chrm	snppos	ref	mat_gtyp	pat_gtyp	c_gtyp	phase	mat_all	pat_all	cA	cC	cG	cT	winning	SymCls	SymPval	BindingSite	cnv
	#NA06986	1	1080920	G	S	W	R	PHASED	G	A	2	0	7	0	M	Sym	0.1796875	1	1.0
	my ($sample_id, $asb_chr,$asb_snppos, $ref, $m_allele, $p_allele, $cA, $cC, $cG, $cT, $winner, $symcls) = (split /\t/)[0,1,2,3,8,9,10,11,12,13,14,15];
	
	next if(!$asb_chr);
	next if(!$asb_snppos);
	if($sym_variants){
		next unless ($symcls =~ /Sym/);
	}else{
		next unless ($symcls =~ /Asym/);
	}
	my $asb_coordinate_id =  $asb_chr . '-' . $asb_snppos;

	my $allele_ancestral = $variants_1kg{$asb_coordinate_id}{ANC};
	my $allele_derived = $variants_1kg{$asb_coordinate_id}{DER}; 
	next unless($allele_ancestral);
	next unless($allele_derived);
#	my $ancestral = $ref;
#	my $derived;
#	if($ref eq $m_allele){
#		$derived = $p_allele;
#	}elsif($ref eq $p_allele){
#		$derived = $m_allele;
#	}else{
#		next;
#	}	
	
	if( ($m_allele eq 'None') or ($p_allele eq 'None')  ){ #none happens for unphased genotypes (imputed)
	}else{
		if( ($allele_ancestral eq $m_allele)  && ($allele_derived eq $p_allele) ){		
		}elsif( ($allele_ancestral eq $p_allele)  && ($allele_derived eq $m_allele)  ){	
		}else{
			print STDERR "$asb_coordinate_id: (ancestral/derived and mat/pat don't match: ($allele_ancestral, $allele_derived) and ($m_allele, $p_allele)\n";
			next;
		}		
	}
	my $read_data = join(",", $allele_ancestral,$allele_derived, $cA, $cC, $cG, $cT);
	$position2sample2readdepth{$asb_coordinate_id}{$sample_id} = $read_data;
}
close $instream;


#now I select variants (create flag recurrent if you need)
#and save, for each VDR-BV that is a QTL, the direction of binding and the direction of expression 
my %data_entries;
open ($instream,  q{<}, $input_funseq) or die("Unable to open $input_funseq : $!");
while(<$instream>){
	chomp;
	next if($_ =~ /^\#/);
	next if($_ eq '');
	my $info_motifbr;my $hit_position;
	my $vdr_bv_recur = 'N'; my $vdr_motifbr = 'N';
	
	my ($funseq_chr, $funseq_pos, $funseq_snp_ids, $funseq_info) = (split /\t/)[0,1,2,7];
	my $chr_b37 = $funseq_chr; 
	$chr_b37 =~ s/chr(.*)/$1/; #alleleseq is b37, this is hg19	
	my $funseq_genomic_coord = $chr_b37 . '-' . $funseq_pos;	
	if(!$position2sample2readdepth{$funseq_genomic_coord}){
		print STDERR "$funseq_genomic_coord - $funseq_snp_ids: no ancestral/derived info available for this VDR-BV. Skipping..\n";
		next;
	}
	next if(!$snp2QTL_info{$funseq_genomic_coord});
	
	#tag recur
	$vdr_bv_recur = 'Y' if($funseq_info =~ /RECUR/);
	#next unless ($funseq_info =~ /RECUR/);
	
	#next if the snp id does not match
	#sometimes you have an id in the form rs115138762;rs1362070 from funseq 
	#split
	my @funseq_snp_ids = split(";",$funseq_snp_ids);
	foreach my $funseq_snp_id (@funseq_snp_ids){
		#$check coordinates and snp id
		unless($snp2QTL_info{$funseq_genomic_coord}{$funseq_snp_id}){
			print STDERR "WARNING: rs id: $funseq_snp_id not found at genomic coordinates: $funseq_genomic_coord in the GEUVADIS data. Different build? Skipping..CHECK\n";
			next;
		}
		
		my @info = split(";", $funseq_info);
		#I want to flag rows for which
		#the VDR-BV/eQTL breaks a VDR motif
		foreach my $item (@info){
			$info_motifbr = $item if( ($item =~ /^MOTIFBR/) && ($item =~ /VDR_JASPAR/));				
		}
	
		#if you have motifbr data, you can remove instances where motif pwm change is discordant with phenotype change.
		#if you do not have motifbr data, you cannot do this, and you will have to look at the final results
		#TODO can you use motif br data for other TFs?
		if($info_motifbr){#===================================================================================================
			$vdr_motifbr = 'Y';
			my ($motif_change_type,$motif_change_byTF) = split("=", $info_motifbr); my @motif_change_byTF = split(",", $motif_change_byTF);
			foreach my $item (@motif_change_byTF){
				my ($TF,$TF_motif,$TF_motif_start,$TF_motif_end,$TF_motif_strand,$TF_snp_position,$TF_der_score,$TF_anc_score) =  split("#", $item);
				next unless( $TF_motif eq $motif_name);
				next if(!$TF_motif_strand);
				next unless( ($TF_motif_strand eq '-') or ($TF_motif_strand eq '+') );
				
				my %lob_samples; my %gob_samples; my %discordant_samples; my $BINDING_DIR;
				$lob_samples{$funseq_genomic_coord} = 0; $gob_samples{$funseq_genomic_coord} = 0; $discordant_samples{$funseq_genomic_coord} = 0;
				#Compute gt_score from $TF_anc_score and $TF_alt_score
				#these should already account for ancestral allele in funseq
				my $GT_SCORE = ( $TF_anc_score + 1 ) / ($TF_der_score + 1); #should these all be  LOB by definition? And therefore GT_SCORE > 1 ?
				#Compute ph score and evaluate concordance with gt score at this position
				#if no concordance, it's meaningless to get the DAF
				foreach my $sample (sort keys %{ $position2sample2readdepth{$funseq_genomic_coord} }){
					my ($sample_values, $sample_fc) = process_read_data($position2sample2readdepth{$funseq_genomic_coord}{$sample}, $sample, $funseq_snp_id);
					if( ($GT_SCORE > 1) && ($sample_fc > 1) ){
						$lob_samples{$funseq_genomic_coord} += 1;
					}elsif( ($GT_SCORE > 1) && ($sample_fc < 1)  ){
						$discordant_samples{$funseq_genomic_coord} += 1;
					}elsif( ($GT_SCORE < 1) && ($sample_fc > 1)  ){
						$discordant_samples{$funseq_genomic_coord} += 1;
					}elsif( ($GT_SCORE < 1) && ($sample_fc < 1)  ){
						print STDERR "$funseq_snp_id: MOTIFBR gt-ph concordance, but GT_SCORE is negative?? Check.\n";
						$gob_samples{$funseq_genomic_coord} += 1;
					}else{
						print STDERR "Should not be here: GT score: $GT_SCORE and sample FC: $sample_fc\n";
					}
				}
				#Now I know, for each samples how many are concordant and how many are discordant. 
				#All are LOB, unless something strange is going on, but some will be false.
				if($discordant_samples{$funseq_genomic_coord} > 0){
					if( ($lob_samples{$funseq_genomic_coord} > 0) || ($gob_samples{$funseq_genomic_coord} > 0)  ){
						print STDERR "Warning - motifbr for VDR-BV $funseq_snp_id shows:\n"; 
						print STDERR "concordance (LOB) for $lob_samples{$funseq_genomic_coord} samples,\n";
						print STDERR "concordance (GOB) for $gob_samples{$funseq_genomic_coord} samples and\n"; 
						print STDERR "discordance for $discordant_samples{$funseq_genomic_coord} samples.\n";
						print STDERR "Skip for now. EVALUATE THIS\n";
						next;
					}			
					print STDERR "Warning - motifbr for VDR-BV $funseq_snp_id only shows gt-ph discordance ($discordant_samples{$funseq_genomic_coord}). Skipping..\n";
					next;
				}
				#now I'm left with concordant, label them
				if(  ($lob_samples{$funseq_genomic_coord} > 0 ) && ($gob_samples{$funseq_genomic_coord} > 0)  ){
					print STDERR "Warning VDR-BV $funseq_snp_id - at this position, $lob_samples{$funseq_genomic_coord} sample are LOB and $gob_samples{$funseq_genomic_coord} are GOB.\n";
					if($lob_samples{$funseq_genomic_coord}/$gob_samples{$funseq_genomic_coord}  >= 4){
						print STDERR "Choosing LOB.\n";
						$BINDING_DIR = $LOB;
					}elsif($gob_samples{$funseq_genomic_coord}/$lob_samples{$funseq_genomic_coord}  >= 4){
						print STDERR "Choosing GOB.\n";
						$BINDING_DIR = $GOB;
					}else{
						print STDERR "Skipping\n";						
						next;
					}
				}elsif( ($lob_samples{$funseq_genomic_coord} > 0 ) && ($gob_samples{$funseq_genomic_coord} == 0) ){
					$BINDING_DIR = $LOB;
				}elsif( ($gob_samples{$funseq_genomic_coord} > 0 ) && ($lob_samples{$funseq_genomic_coord} == 0)  ){
					$BINDING_DIR = $GOB;
				}else{
					print STDERR "ERROR VDR-BV $funseq_snp_id $lob_samples{$funseq_genomic_coord} lob samples and  $gob_samples{$funseq_genomic_coord} gob samples - should not be here.\n";
					next;
				}
				
	
				#my %expression_values;
				foreach my $gene_id (keys %{ $snp2QTL_info{$funseq_genomic_coord}{$funseq_snp_id} }){
					my $EFFECT; 
					my $rval = $snp2QTL_info{$funseq_genomic_coord}{$funseq_snp_id}{$gene_id}{RVAL};
					my $pval = $snp2QTL_info{$funseq_genomic_coord}{$funseq_snp_id}{$gene_id}{PVAL};
					my $EXPR_DIR = get_expression_direction($funseq_genomic_coord, $rval);
					if(  ( ($BINDING_DIR eq $LOB) && ($EXPR_DIR eq $LOE) ) || ( ($BINDING_DIR eq $GOB) && ($EXPR_DIR eq $GOE) ) ){
						$EFFECT = $ACTIVATE;
					}elsif(  ( ($BINDING_DIR eq $LOB) && ($EXPR_DIR eq $GOE) ) || ( ($BINDING_DIR eq $GOB) && ($EXPR_DIR eq $LOE) ) ){
						$EFFECT = $REPRESS;		
					}else{
						print STDERR "$funseq_snp_id: directionality problem. Binding dir: $BINDING_DIR; Expression dir: $EXPR_DIR. Skipping..\n";
						next;
					} 	
					my $field_funseq_anno = process_funseq_annotation($funseq_info);
					my $data_line = $funseq_chr    . "\t" . 
									$funseq_pos    . "\t" . 
									$funseq_snp_id . "\t" .
									$vdr_motifbr   . "\t" .
									$vdr_bv_recur  . "\t" .
									$gene_id       . "\t" . 
									$rval          . "\t" .
									$pval          . "\t" . 
									$BINDING_DIR   . "\t" . 
									$EXPR_DIR      . "\t" . 
									$EFFECT        . "\t" .
									$field_funseq_anno;
					#print $data_line, "\n";
					$data_entries{$data_line} = 1;
				}
			}	
		}else{ #no motifbr - the variant can be anywhere
			#get lob and gob label for all samples
			my %lob_samples; my %gob_samples; my $BINDING_DIR;
			$lob_samples{$funseq_genomic_coord} = 0; $gob_samples{$funseq_genomic_coord} = 0;
			foreach my $sample (sort keys %{ $position2sample2readdepth{$funseq_genomic_coord} }){
				my ($values, $fc) = process_read_data($position2sample2readdepth{$funseq_genomic_coord}{$sample}, $sample, $funseq_snp_id);
				if($fc > 1){
					$lob_samples{$funseq_genomic_coord} += 1;
				}elsif($fc < 1){
					$gob_samples{$funseq_genomic_coord} += 1;
				}elsif($fc == 1){
					print STDERR "Warning: VDR-BV $funseq_snp_id, sample $sample - phenotype fc is $fc. Skipping sample..\n";
					next;
				}else{
					print "Warning fold change value $fc is empty or undefined. Skipping..\n";
					next;
				}
			}	
			#decide whether to keep or discard
			if(  ($lob_samples{$funseq_genomic_coord} > 0 ) && ($gob_samples{$funseq_genomic_coord} > 0)  ){
				print STDERR "Warning VDR-BV $funseq_snp_id - at this position, $lob_samples{$funseq_genomic_coord} sample are LOB and $gob_samples{$funseq_genomic_coord} are GOB.\n";
				if($lob_samples{$funseq_genomic_coord}/$gob_samples{$funseq_genomic_coord}  >= 4){
					print STDERR "Choosing LOB.\n";
					$BINDING_DIR = $LOB;
				}elsif($gob_samples{$funseq_genomic_coord}/$lob_samples{$funseq_genomic_coord}  >= 4){
					print STDERR "Choosing GOB.\n";
					$BINDING_DIR = $GOB;
				}else{
					print STDERR "Skipping\n";						
					next;
				}
			}elsif( ($lob_samples{$funseq_genomic_coord} > 0 ) && ($gob_samples{$funseq_genomic_coord} == 0) ){
				$BINDING_DIR = $LOB;
			}elsif( ($gob_samples{$funseq_genomic_coord} > 0 ) && ($lob_samples{$funseq_genomic_coord} == 0)  ){
				$BINDING_DIR = $GOB;
			}else{
				print STDERR "ERROR VDR-BV $funseq_snp_id $lob_samples{$funseq_genomic_coord} lob samples and  $gob_samples{$funseq_genomic_coord} gob samples - should not be here.\n";
				next;
			}	

			#my %expression_values;
			foreach my $gene_id (keys %{ $snp2QTL_info{$funseq_genomic_coord}{$funseq_snp_id} }){
				my $EFFECT; 
				my $rval = $snp2QTL_info{$funseq_genomic_coord}{$funseq_snp_id}{$gene_id}{RVAL};
				my $pval = $snp2QTL_info{$funseq_genomic_coord}{$funseq_snp_id}{$gene_id}{PVAL};
				my $EXPR_DIR = get_expression_direction($funseq_genomic_coord, $rval);
				if(  ( ($BINDING_DIR eq $LOB) && ($EXPR_DIR eq $LOE) ) || ( ($BINDING_DIR eq $GOB) && ($EXPR_DIR eq $GOE) ) ){
					$EFFECT = $ACTIVATE;
				}elsif(  ( ($BINDING_DIR eq $LOB) && ($EXPR_DIR eq $GOE) ) || ( ($BINDING_DIR eq $GOB) && ($EXPR_DIR eq $LOE) ) ){
					$EFFECT = $REPRESS;		
				}else{
						print STDERR "$funseq_snp_id: directionality problem. Binding dir: $BINDING_DIR; Expression dir: $EXPR_DIR. Skipping..\n";
					next;
				} 	
				my $field_funseq_anno = process_funseq_annotation($funseq_info);
				my $data_line = $funseq_chr    . "\t" . 
								$funseq_pos    . "\t" . 
								$funseq_snp_id . "\t" .
								$vdr_motifbr   . "\t" .
								$vdr_bv_recur  . "\t" .
								$gene_id       . "\t" . 
								$rval          . "\t" .
								$pval          . "\t" . 
								$BINDING_DIR   . "\t" . 
								$EXPR_DIR      . "\t" . 
								$EFFECT        . "\t" .
								$field_funseq_anno;
					#print $data_line, "\n";
					$data_entries{$data_line} = 1;
			}			
		}	
	}
}
close $instream;

open (my $outstream,  q{>}, $output) or die("Unable to open $output : $!");
print $outstream "VDRBV_CHR\tVDRBV_POS\tVDRBV_ID\tisVDRRXR-MOTIFBR\tisVDR-hcBV\teQTL_GENE_ID\teQTL_RVAL\teQTL_PVAL\tVDR_BINDING_DIR\tGENE_EXPR_DIR\tEFFECT\tFUNSEQ_GENE\tFUNSEQ_NC_ANNOT\tFUNSEQ_HOT\tFUNSEQ_TFP\tFUNSEQ_TFM\tFUNSEQ_PROTNET_HUB\tFUNSEQ_GERP\tFUNSEQ_RECUR\n";
foreach my $item (sort keys %data_entries){
	print $outstream $item, "\n";
}
close $outstream;

#============================================
sub process_read_data{
	my ($read_data_string, $lcl_sample, $id) = @_;
	
	my ($ancestral_ref,$ancestral_alt,$cA,$cC,$cG,$cT) = split(",", $read_data_string);
	my $anc_count; my $der_count;
	my $der_base = $ancestral_alt;
	my $anc_base = $ancestral_ref;
	
	if( ($cA eq 0) && ($cC eq 0) && ($cG eq 0) && ($cT eq 0) ){
		print STDERR "process_read_data(): error asb snp: $id. (LCL sample: $lcl_sample) - 0 read coverage for all possible nucleotides at this event ($cA,$cC,$cG,$cT). Aborting..\n";
		exit -1;
	}

	#get reference read coverage
	if($anc_base eq 'A'){
		$anc_count = $cA;
	}elsif($anc_base eq 'C'){
		$anc_count = $cC;
	}elsif($anc_base eq 'G'){
		$anc_count = $cG;
	}elsif($anc_base eq 'T'){
		$anc_count = $cT;
	}else{
		print STDERR "process_read_data(): error asb snp: $id. (LCL sample: $lcl_sample) - reference base not recognised: $anc_base. Aborting..\n";
		exit -1;
	}
	
	if($der_base eq 'A'){
		$der_count = $cA;
	}elsif($der_base eq 'C'){
		$der_count = $cC;
	}elsif($der_base eq 'G'){
		$der_count = $cG;
	}elsif($der_base eq 'T'){
		$der_count = $cT;
	}else{
		print STDERR "process_read_data(): error asb snp: $id. (LCL sample: $lcl_sample) - alternate base not recognised: $der_base. Aborting..\n";
		exit -1;
	}
	my $fold_change = ( $anc_count + 1 ) / ( $der_count + 1 );
	
	#I will let R carry out the division.
	#I will save the value (REF)-(ALT)
	my $value_string = $anc_count . '-' . $der_count;
	return ($value_string, $fold_change);
}


#query the hash for the ancestral and derived allele
#then assign
#loss of expression, gain of expression
#$LOE = 'loe'; #ancestral -> derived reduces expression  
#$GOE = 'goe'; #ancestral -> derived increases expression

sub get_expression_direction{
	my ($coord_id, my $rvalue) = @_;
	#The sign denotes the direction of the nonreference allele (i.e. rvalue<0 means that nonreference (alternate) allele has lower expression)
	#If not swap alternate = derived; if minus, LOE, if not minus GOE 
	#If swap, alternate = ancestral; if minus, GOE; if not minus, LOE

#	if($variants_1kg{$coord_id}{SWAP}){ #the non reference (alternate) is the ancestral TEMP
#		if($rvalue < 0){   
#			return $GOE;
#		}else{
#			return $LOE;
#		}
#	}else{ #the non reference (alternate) is the derived; no swap
#		if($rvalue < 0){
#			return $LOE;
#		}else{
#			return $GOE;
#		}
#	}	
	
			if($rvalue < 0){   
			return $GOE;
		}else{
			return $LOE;
		}
}

#I have something like
#SAMPLE=interestingHets_NA19213_EBLfiltered_hg19;GERP=-3.66;CDS=No;HUB=FAM200B:REG(0.880),FBXL5:PPI(0.361);NCENC=DHS(MCV-42|chr4:15755985-15756135),Enhancer(chmm/segway|chr4:15755162-15756284),Enhancer(drm|chr4:15756000-15756200),TFM(VDR|VDR_JASPAR|chr4:15756098-15756113),TFM(VDR|VDR_XXmotif|chr4:15756098-15756113),TFM(VDR|VDR_dreme|chr4:15756098-15756105),TFP(FOS|chr4:15755301-15756613),TFP(GATA2|chr4:15755241-15756612),TFP(MAX|chr4:15755846-15756384),TFP(MYC|chr4:15755856-15756339),TFP(VDR|chr4:15755910-15756250);MOTIFBR=VDR#VDR_dreme#15756098#15756105#-#4#0.000000#0.476400,VDR#VDR_JASPAR#15756098#15756113#-#12#0.000#1.000,VDR#VDR_XXmotif#15756098#15756113#-#12#0.02457#0.73562;GENE=FAM200B(Distal)[pearson(H3K27ac):0.918602,pearson(H3K4me1):0.912105],FBXL5(Distal)[pearson(H3K27ac):0.875722,pearson(H3K4me1):0.873539];NCDS=1.67156771704517

#or

#SAMPLE=interestingHets_NA19248_EBLfiltered_hg19;GERP=0.375;CDS=No;NCENC=DHS(MCV-31|chr7:7298040-7298190),Enhancer(chmm/segway|chr7:7296800-7298400),TFM(VDR|VDR_JASPAR|chr7:7298122-7298137),TFP(BATF|chr7:7297738-7298801),TFP(BCL11A|chr7:7297625-7298314),TFP(IRF4|chr7:7297670-7298802),TFP(MAX|chr7:7297536-7298947),TFP(MEF2A|chr7:7297922-7298725),TFP(MYC|chr7:7297025-7298988),TFP(NFKB1|chr7:7297551-7298769),TFP(PAX5|chr7:7297634-7298780),TFP(PAX5|chr7:7297644-7298634),TFP(PAX5|chr7:7297659-7298342),TFP(PAX5|chr7:7297918-7298625),TFP(POU2F2|chr7:7297992-7298671),TFP(RXRA|chr7:7297871-7298750),TFP(SPI1|chr7:7297558-7298322),TFP(SPI1|chr7:7297604-7298306),TFP(VDR|chr7:7298055-7298256);HOT=Gm12878;MOTIFBR=VDR#VDR_JASPAR#7298122#7298137#+#12#0.000#1.000;NCDS=1.788999949

#SAMPLE=interestingHets_NA19213_EBLfiltered_hg19;GERP=-3.66;CDS=No;HUB=FAM200B:REG(0.880),FBXL5:PPI(0.361);NCENC=DHS(MCV-42|chr4:15755985-15756135),Enhancer(chmm/segway|chr4:15755162-15756284),Enhancer(drm|chr4:15756000-15756200),TFM(VDR|VDR_JASPAR|chr4:15756098-15756113),TFM(VDR|VDR_XXmotif|chr4:15756098-15756113),TFM(VDR|VDR_dreme|chr4:15756098-15756105),TFP(FOS|chr4:15755301-15756613),TFP(GATA2|chr4:15755241-15756612),TFP(MAX|chr4:15755846-15756384),TFP(MYC|chr4:15755856-15756339),TFP(VDR|chr4:15755910-15756250);MOTIFBR=VDR#VDR_dreme#15756098#15756105#-#4#0.000000#0.476400,VDR#VDR_JASPAR#15756098#15756113#-#12#0.000#1.000,VDR#VDR_XXmotif#15756098#15756113#-#12#0.02457#0.73562;GENE=FAM200B(Distal)[pearson(H3K27ac):0.918602,pearson(H3K4me1):0.912105],FBXL5(Distal)[pearson(H3K27ac):0.875722,pearson(H3K4me1):0.873539];NCDS=1.67156771704517

#I want
#FUNSEQ_GENE\tFUNSEQ_NC_ANNOT\tFUNSEQ_HOT\tFUNSEQ_TFP\tFUNSEQ_TFM\tFUNSEQ_PROTNET_HUB\tFUNSEQ_GERP\tFUNSEQ_RECUR\n

sub process_funseq_annotation{
	my ($funseq_string) = @_;
	#output fields
	my $info_gerp; my $info_hub; my $info_ncenc; my $info_hot; my $info_gene; my $info_recur;my $info_tfp; my $info_tfm;
	my $out_gerp; my $out_hub; my $out_ncenc; my $out_hot; my $out_gene; my $out_recur; my $out_tfp = ''; my $out_tfm = '';
	my @info = split(";",$funseq_string);
	
	foreach my $item (@info){
		$info_gerp = $item if( ($item =~ /^GERP/) );
		$info_hub = $item if( ($item =~ /^HUB/) ) ;				
		$info_ncenc = $item if( ($item =~ /^NCENC/) );
		$info_hot = $item if( ($item =~ /^HOT/) );
		$info_gene = $item if( ($item =~ /^GENE/) );
		$info_recur = $item if( ($item =~ /^RECUR/) );		
	}
	#gerp-----------------------------------------------------
	if($info_gerp =~ /GERP=(.*)/){
		$out_gerp = $1;
	}else{
		$out_gerp = '';
	}
	#hub-----------------------------------------------------
	if($info_hub){
		my %entries; my @entries;
		my ($header, $data) = split('=',$info_hub);
		my @info_hub = split(',', $data);
		foreach my $item (@info_hub){
			if($item =~ /(.+)\:(.+)/){
				$entries{$1} = 1;
			}else{
				next;
			}				
		}
		foreach my $item (sort keys %entries){
			push(@entries, $item);	
		}
		$out_hub = join(',',@entries);
	}else{
		$out_hub = '';
	}
	#ncenc-----------------------------------------------------
#qua devi fare attenzione
#NCENC=
#DHS(MCV-42|chr4:15755985-15756135),
#Enhancer(chmm/segway|chr4:15755162-15756284),
#Enhancer(drm|chr4:15756000-15756200),
#TFM(VDR|VDR_JASPAR|chr4:15756098-15756113),
#TFM(VDR|VDR_XXmotif|chr4:15756098-15756113),
#TFM(VDR|VDR_dreme|chr4:15756098-15756105),
#TFP(FOS|chr4:15755301-15756613),
#TFP(GATA2|chr4:15755241-15756612),
#TFP(MAX|chr4:15755846-15756384),
#TFP(MYC|chr4:15755856-15756339),
#TFP(VDR|chr4:15755910-15756250);	
	
	if($info_ncenc){
		my %entries; my @entries;
		my %tfm_entries; my @tfm_entries;
		my %tfp_entries; my @tfp_entries;		
		
		my ($header, $data) = split('=',$info_ncenc);
		my @info_ncenc = split(',', $data);
		
		#if there is encode peak/motif annotation, process and store separately
		foreach my $item (@info_ncenc){
			if( ($item =~ /^TFM\((.*)\|chr(.*)\)/) ){#TFM(VDR|VDR_JASPAR|chr4:15756098-15756113)
				$tfm_entries{$1} = 1;
			}
			if( ($item =~ /^TFP\((.*)\|chr(.*)\)/) ){#TFP(GATA2|chr4:15755241-15756612)
				$tfp_entries{$1} = 1;				
			}
		}
		foreach my $item (sort keys %tfm_entries){
			push(@tfm_entries, $item);	
		}		
	 	$out_tfm = join(',',@tfm_entries);	
	 	
		foreach my $item (sort keys %tfp_entries){
			push(@tfp_entries, $item);	
		}		
	 	$out_tfp = join(',',@tfp_entries);	
		
		#back to ncenc
		foreach my $item (@info_ncenc){
			if($item =~ /(.+)\(.+\)/){
				$entries{$1} = 1;
			}else{
				next;
			}				
		}
		foreach my $item (sort keys %entries){
			push(@entries, $item);	
		}		
	 	$out_ncenc = join(',',@entries);
	}else{
		$out_ncenc = '';
	}

	#hot-----------------------------------------------------
	if($info_hot){
		$out_hot = 'Y';	
	}else{
		$out_hot = '';
	}
	#GENE=FAM200B(Distal)[pearson(H3K27ac):0.918602,pearson(H3K4me1):0.912105],FBXL5(Distal)[pearson(H3K27ac):0.875722,pearson(H3K4me1):0.873539]
	#gene-----------------------------------------------------
	if($info_gene){
		my %entries; my @entries;
		my ($header, $data) = split('=',$info_gene);
		my @info_gene = split('],', $data);
		
		foreach my $item (@info_gene){
			if($item =~ /(.+)\((.+)\)\[.*/){
				my $gene = $1 . '(' . $2 . ')';
				$entries{$gene} = 1;
			}elsif($item =~ /(.+)\((.+)\)/){
				my $gene = $1 . '(' . $2 . ')';
				$entries{$gene} = 1;
			}else{
				next;
			}				
		}
		foreach my $item (sort keys %entries){
			push(@entries, $item);	
		}		
	 	$out_gene = join(',',@entries);		
	}else{
		$out_gene = '';
	}
	
	#recur-----------------------------------------------------
	if($info_recur){
		$out_recur = 'Y';
	}else{
		$out_recur = '';
	}
	my $out_string = $out_gene  . "\t" .
					 $out_ncenc . "\t" .
					 $out_hot   . "\t" .
					 $out_tfp   . "\t" .
					 $out_tfm   . "\t" .					  
					 $out_hub   . "\t" .
					 $out_gerp  . "\t" .
					 $out_recur;

	return $out_string;
}