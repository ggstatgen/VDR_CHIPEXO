#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use List::Util qw(min max);


#run this after running impute2. This is to filter impute genotypes with low reliability
#USELESS - you should be able to do everything with qctool;
#http://www.well.ox.ac.uk/~gav/qctool/#tutorial


#on the IMPUTE mailing list I found the following:
#From: Oxford Statistical Genetics Software [mailto:OXSTATGEN@JISCMAIL.AC.UK] On Behalf Of Vince Forgetta
#Sent: Tuesday, October 01, 2013 11:32 AM
#To: OXSTATGEN@JISCMAIL.AC.UK
#Subject: Re: [OXSTATGEN] Post-imputation QC
#Hi Majid,
#
#As Sarah says, INFO is the primary filter we use to filter for quality after imputation.
#
#
#Typically we discard SNPs with INFO < 0.4 *and* HWE p < 1e-6. However, these are both subject to change depending on aims of the study and other factors. I am by no means an #expert and defer to others here for additional and wiser comments. But from my experience, our group is considering a higher INFO threshold (e.g. INFO < 0.8) for rare variants #(MAF < 1%) to ensure we get good quality imputation for these low frequency variants.
#
#
#Hope this helps,
#
#
#
#Vince 

#2
#I agree with all the comments below. There is no single correct set of filters and you need to try adjusting the filters as needed for you particular study. So start with no #filters and gradually apply them.

#Jonathan

#3 http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3410659/

#For now I will remove all snps with imputed quality < 0.4