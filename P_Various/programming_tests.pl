#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;


my $infile = "file.txt";

open (my $instream,  q{<}, $infile) or die("Unable to open $infile : $!");
close $infile;
