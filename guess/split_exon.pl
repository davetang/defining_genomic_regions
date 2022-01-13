#!/usr/bin/env perl

use strict;
use warnings;

while(<>){
   chomp;
   next if /^chrom/;
   my ($chr, $starts, $ends) = split(/\t/);
   my @starts = split(/,/, $starts);
   my @ends = split(/,/, $ends);
   foreach my $i (0..$#starts){
      print join("\t", $chr, $starts[$i], $ends[$i]), "\n";
   }
}

exit();

