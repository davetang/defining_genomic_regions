#!/bin/env perl

use strict;
use warnings;

my $v = 19;

my $exon_file       = "gencode_v${v}_exon_merged.bed.gz";
my $intergenic_file = "gencode_v${v}_intergenic.bed.gz";
my $intron_file     = "gencode_v${v}_intron.bed.gz";

my $exon_coverage       = coverage($exon_file);
my $intergenic_coverage = coverage($intergenic_file);
my $intron_coverage     = coverage($intron_file);

my $total = $exon_coverage + $intergenic_coverage + $intron_coverage;

printf "Exon: %.2f\n", $exon_coverage*100/$total;
printf "Intron: %.2f\n", $intron_coverage*100/$total;
printf "Intergenic: %.2f\n", $intergenic_coverage*100/$total;

sub coverage {
   my ($infile) = @_;
   my $coverage = 0;
   open(IN,'-|',"zcat $infile") || die "Could not open $infile: $!\n";
   while(<IN>){
      chomp;
      my ($chr, $start, $end) = split(/\t/);
      my $c = $end - $start;
      $coverage += $c;
   }
   close(IN);
   return($coverage);
}

exit(0);
