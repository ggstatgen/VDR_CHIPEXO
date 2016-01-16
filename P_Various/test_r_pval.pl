#!/usr/bin/perl
use strict;
use warnings;
use Statistics::R;

my @pvals = (0.266845129, 0.353675275,0.030585143,0.791285527,0.60805832,0.433466716,0.039369439,0.48918122,0.740626586);
my $R = Statistics::R->new();
$R->set('pvals', \@pvals);
$R->run(q`adjustp <- p.adjust(pvals, method = "BH", n = length(pvals))`);
my $corrected_pvals = $R->get('adjustp');
$R->stop();

my $counter = 0;
print "PVAL\tCORR_PVAL\n";
foreach my $item (@pvals){
	print $item . "\t" . $$corrected_pvals[$counter] . "\n";
	$counter++; 
} 