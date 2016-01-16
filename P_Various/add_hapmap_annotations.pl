#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#14/6/2013
#I used the GATK functions to convert my hapmap files into vcfs.
#However these vcfs don't have any annotation
#hapmap also provides allele frequencies divided by population
#I use this script to enrich one .vcf file obtained from hapmap with the frequencies
#Frequencies should be the same release (using r27 atm)

