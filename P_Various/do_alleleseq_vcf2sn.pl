#!/usr/bin/perl -w
use strict;
use warnings;
#use fralib;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

###
#MODIFICATO GIU 30/1/14
###
#THIS EXPECTS A VCF WITH THREE GENOTYPE COLUMNS: MOTHER	FATHER	CHILD
#$Id: vcf_2_snp.pl,v 1.3 2013/01/19 19:34:00 yk336 Exp $
#
# convert vcf format to file needed by alleleseq pipeline
#


=head1 NAME

vcf2snp

=head1 SYNOPSIS

 vcf2snp <vcf> > out.txt

  -h help
  
  Convert vcf format to file needed by alleleseq pipeline
  GIUSEPPE: it expects a vcf containing genotype for: MOTHER,FATHER,CHILD
                
 
  Example:
     vcf2snp snv.vcf > out.txt
  
=head1 DESCRIPTION

=cut

#option variables
my $help;

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help) || scalar(@ARGV)!=1){
    if ($help){
        pod2usage(-verbose => 2);
    }
    else{
        pod2usage(1);
    }
}
my $fn = shift;

my %iupac2code = (
                  'A' => 1,
                  'C' => 2,
                  'G' => 4,
                  'T' => 8,
                  'R' => 1|4,
                  'Y' => 2|8,
                  'S' => 2|4,
                  'W' => 1|8,
                  'K' => 4|8,
                  'M' => 1|2,
                  'B' => 2|4|8,
                  'D' => 1|4|8,
                  'H' => 1|2|8,
                  'V' => 1|2|4,
                  'N' => 1|2|4|8,
                  );

my %code2iupac = (
                  1 => 'A',
                  2 => 'C',
                  4 => 'G',
                  8 => 'T',
                  1|4 => 'R',
                  2|8 => 'Y',
                  2|4 => 'S',
                  1|8 => 'W',
                  4|8 => 'K',
                  1|2 => 'M',
                  2|4|8 => 'B',
                  1|4|8 => 'D',
                  1|2|8 => 'H',
                  1|2|4 => 'V',
                  1|2|4|8 => 'N',
    );


my $fh; 
open($fh, "<$fn") || die "cannot open the file $fn!";
while (my $l = <$fh>){
    chomp($l);
    next if ($l =~ /^\#/);

    my @t = split(/\t/, $l);
    my $chr = $t[0]; 
    my $pos = $t[1];
    my $ref = $t[3];
    my $alt = $t[4];

#    next if ($t[6] ne "PASS"); #GIUSEPPE
#    next if (length($ref) != length($alt));
    next if (length($ref) > 1);
    my @alts = split(/,/, $alt);
    next if (length($alts[0]) > 1); # child

    ### now @alts 1st is ref
    unshift(@alts, $ref);

    my $format = $t[8];
    #GIUSEPPE------------------------
    #my $chd = $t[9];
    #my $pat = $t[10];
    #my $mat = $t[11];
    my $mat = $t[9];
    my $pat = $t[10];
    my $chd = $t[11];
    #GIUSEPPE-----------------------
    my @fmts = split(/:/, $format);
    my $idx = 0;
    for (@fmts) {
		last if ($_ eq "GT");
		$idx++;
    }

    my ($c1, $c2, $phased) = &gt($chd, $idx);
    my ($p1, $p2) = &gt($pat, $idx);
    my ($m1, $m2) = &gt($mat, $idx);

    next if ($c1 eq "." || $c2 eq "."); # child unknown

    # mother first?
    my $ca = $alts[$c2] . $alts[$c1];
    my $pa = ($p2 eq "." ? "." : $alts[$p2]) . ($p1 eq "." ? "." : $alts[$p1]);
    my $ma = ($m2 eq "." ? "." : $alts[$m2]) . ($m1 eq "." ? "." : $alts[$m1]);

    $pa = &guess($ca, $ma) if ($pa =~ /\./);
    $ma = &guess($ca, $pa) if ($ma =~ /\./);
    
    if ($pa =~ /\./ || $ma =~ /\./) {
	print STDERR "$l\n";
	next;
    }

    my $mutant = &is_mutant($ca, $pa, $ma);

    my $status;
    if ($mutant){
		$status = "MUTANT";
    }elsif ($c1 == $c2) {
		$status = "HOMO";
    }elsif ($phased) {
		$status = "PHASED";
    }elsif ($c1 != $c2) {
		$status = "HETERO";
    }else{
		$status = "UNKNOWN";
    }

=pod
    print join("\t",
	       $chr,
	       $pos,
	       $ref,
	       $alt,
	       $chd,
	       $pat,
	       $mat,
	       @alts,
	       $idx,
	       $c1, $c2,
	       $p1, $p2,
	       $m1, $m2,
	       $ca, $pa, $ma,
	       ">$mutant<",
	       $status,
	       ), "\n";
=cut

#=pod
    #GIUSEPPE COMMENTED THE FOLLOWING LINE	
    $chr =~ s/^chr//;
    print join("\t",
	       $chr,
	       $pos,
	       $ref,
	       $ma,
	       $pa,
	       $ca,
	       $status,
	       ), "\n";
#=cut

}


sub gt {
    my ($s, $idx, ) = @_;

    my @t = split(/:/, $s);
    my $gt = $t[$idx];

#    return (-1, -1) unless ($gt =~ /\|/); # not phased

    my @gts = split(/[\|\/]/, $gt);

    my $phased = 0;
    $phased = 1 if ($gt =~ /\|/);

    return ($gts[0], $gts[1], $phased);
}

### use child and one parent to guess the other parent
sub guess {
    my ($ca, $pa, ) = @_;

    my $r;

    return $ca if ($pa =~ /\./);

    my ($c1, $c2) = split('', $ca);
    my ($p1, $p2) = split('', $pa);

    my $o1 = $iupac2code{$c1} ^ $iupac2code{$p1};
    my $o2 = $iupac2code{$c2} ^ $iupac2code{$p2};

    if ($o1 || $o2) {
	if ($o1) {
	    $r = $c1 . $p1;
	} else {
	    $r = $c2 . $p2;
	}
    } else {
	$r = $ca;
    }
    return $r;
}


sub is_mutant {
    my ($ca, $pa, $ma, ) = @_;

    my ($c1, $c2) = split('', $ca);
    my ($p1, $p2) = split('', $pa);
    my ($m1, $m2) = split('', $ma);

    my $c = $code2iupac{ $iupac2code{$c1} | $iupac2code{$c2} };
    my $p = $code2iupac{ $iupac2code{$p1} | $iupac2code{$p2} };
    my $m = $code2iupac{ $iupac2code{$m1} | $iupac2code{$m2} };

    return &mutant($p, $m, $c);
}

#
# for phasing
#

sub homo {
    my $ltr = shift;

    my $n = $iupac2code{$ltr};
    return $n==1 || $n==2 || $n==4 || $n==8;  
    # note: "or" does not work here
    # due to low precedence than "="
}

sub hetero {
    my $ltr = shift;

    return !homo($ltr);
}

sub comp {
    my ($a, $b) = @_;

    return $iupac2code{$a}^$iupac2code{$b};
}

sub mutant {
    my ($f, $m, $c) = @_;
    my ($m1, $m2);

    my ($fc, $mc, $cc);

    # check if child has any alleles not present in the parents, 
    # indicating
    # a mutation or sequencing error'''
    
    $fc = $iupac2code{$f};
    $mc = $iupac2code{$m};
    $cc = $iupac2code{$c};

    $m1 = $cc & ~($fc | $mc); # outside fc and mc

    $m2 = !( ($cc & $fc) && ($cc & $mc) );
   
#    print "==>>$m1<\t>$m2<\n";

    return ($m1 || $m2);
}

sub a_and_b {
    my ($a, $b) = @_;

    return $iupac2code{$a} & $iupac2code{$b};
}
