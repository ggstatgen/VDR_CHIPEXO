#!/usr/bin/perl -w
# Runs the liftover tool on a VCF and properly handles the output
#19/6/13 obtained from the GATK site, modified by Giuseppe
#does validation before lifting over
use strict;
use warnings;
use Getopt::Long;
use File::Basename;#giu

my $in = undef;
my $gatk = undef;
my $chain = undef;
my $newRef = undef;
my $oldRef = undef;
my $out = undef;
my $tmp = "/tmp";
my $recordOriginalLocation = 0;
GetOptions( "vcf=s" => \$in,
	    "gatk=s" => \$gatk,
	    "chain=s" => \$chain,
	    "newRef=s" => \$newRef,
	    "oldRef=s" => \$oldRef,
            "out=s" => \$out,
	    "tmp=s" => \$tmp,
	    "recordOriginalLocation" => \$recordOriginalLocation);

if ( !$in || !$gatk || !$chain || !$newRef || !$oldRef || !$out ) {
    print "Usage: liftOverVCF.pl\n\t-vcf \t\t<input vcf>\n\t-gatk \t\t<path to gatk trunk>\n\t-chain \t\t<chain file>\n\t-newRef \t<new reference (.fa or .fasta)>\n\t-oldRef \t<old reference (.fa or .fasta)>\n\t-out \t\t<output vcf>\n\t-tmp \t\t<temp file location; defaults to /tmp>\n\t-recordOriginalLocation \t\t<Should we record what the original location was in the INFO field?; defaults to false>\n";
    print "Example: ./liftOverVCF.pl\n\t-vcf /humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/1kg_snp_validation/all_validation_batches.b36.vcf\n\t-chain b36ToHg19.broad.over.chain\n\t-out lifted.hg19.vcf\n\t-gatk /humgen/gsa-scr1/ebanks/Sting_dev\n\t-newRef /seq/references/Homo_sapiens_assembly19/v0/hg19.fa\n\t-oldRef /humgen/1kg/reference/human_b36_both.fa\n";
    exit(1);
}
#giu - get basename for $newRef
my($bn_newref, $dir_nr) = fileparse($newRef);
$bn_newref =~ s/(.*)\..*/$1/;

my $PKEY="/net/isi-scratch/giuseppe/GATK_RESOURCES/giuseppe.gallone_dpag.ox.ac.uk.key";
my $DBSNP_hg18 = "/net/isi-scratch/giuseppe/GATK_RESOURCES/hg18/dbsnp_137.hg18.vcf";

# generate a random number
my $random_number = rand();
my $tmp_prefix = "$tmp/$random_number";
print "Writing temporary files to prefix: $tmp_prefix\n";
my $unsorted_vcf = "$tmp_prefix.unsorted.vcf";

#validate the vcf file (Giuseppe)
#print "Validating the vcf...";
# removed --dbsnp $DBSNP_hg18
#my $validation_out = "java -Xmx2g -jar $gatk/GenomeAnalysisTK.jar -T ValidateVariants -R $oldRef --variant $in --warnOnErrors -et NO_ET -K $PKEY";
#system($validation_out) == 0 or quit("The validation step failed.  Please correct the necessary errors before retrying.");

# lift over the file
print "Lifting over the vcf...";
my $cmd = "java -Xmx2g -jar $gatk/GenomeAnalysisTK.jar -T LiftoverVariants -R $oldRef -V:variant $in -o $unsorted_vcf -chain $chain -dict $dir_nr$bn_newref.dict -et NO_ET -K $PKEY --logging_level INFO";
if ($recordOriginalLocation) {
  $cmd .= " -recordOriginalLocation";
}
system($cmd) == 0 or quit("The liftover step failed.  Please correct the necessary errors before retrying.");

# we need to sort the lifted over file now
print "\nRe-sorting the vcf...\n";
my $sorted_vcf = "$tmp_prefix.sorted.vcf";
open(SORTED, ">$sorted_vcf") or die "can't open $sorted_vcf: $!";

# write the header
open(UNSORTED, "< $unsorted_vcf") or die "can't open $unsorted_vcf: $!";
my $inHeader = 1;
while ( $inHeader == 1 ){
	my $line = <UNSORTED>;
    if ( $line !~ m/^#/ ) {
		$inHeader = 0;
    }else{
		print SORTED "$line";
    }
}
close(UNSORTED);
close(SORTED);

$cmd = "grep \"^#\" -v $unsorted_vcf | sort -n -k2 -T $tmp | $gatk/perl/sortByRef.pl --tmp $tmp - $newRef.fai >> $sorted_vcf";
system($cmd) == 0 or quit("The sorting step failed.  Please correct the necessary errors before retrying.");

# Filter the VCF for bad records
print "\nFixing/removing bad records...\n";
$cmd = "java -jar $gatk/GenomeAnalysisTK.jar -T FilterLiftedVariants -R $newRef -V:variant $sorted_vcf -o $out -et NO_ET -K $PKEY";
system($cmd) == 0 or quit("The filtering step failed.  Please correct the necessary errors before retrying.");

# clean up
unlink $unsorted_vcf;
unlink $sorted_vcf;
my $sorted_index = "$sorted_vcf.idx";
unlink $sorted_index;

print "\nDone!\n";

sub quit {
    print "\n$_[0]\n";
    exit(1);
}
