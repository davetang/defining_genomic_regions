#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('h:s:b:', \%opts);

if ($opts{'h'} ||
    !exists $opts{'s'} ||
    !exists $opts{'b'}
){
   usage();
}

my $genome = $opts{'s'};
my %sizes = store_size($genome);
my $bed = $opts{'b'};

my $fh;
if ($bed =~ /\.gz$/){
   open($fh, '-|', "gunzip -c $bed") or die "Could not open $bed $!\n";
} else {
   open($fh, '<', $bed) or die "Could not open $bed $!\n";
}

while(<$fh>){
   chomp;
   next if /^browser/ || /^track/;
   my ($chr, $start, $end, @rest) = split(/\t/);
   if (exists $sizes{$chr}){
      my $size = $sizes{$chr};
      if (++$end > $size){
         warn("$chr:$start-$end greater than $size\n");
      }
   } else {
      die("$chr does not exist in $genome\n");
   }
}
close($fh);

sub usage {
print STDERR <<EOF;
Usage: $0 -s hg19.genome.gz -b my.bed

Where:   -s         chromosome sizes
         -b         BED file
         -h         this helpful usage message

EOF
exit();
}

sub store_size {
   my ($infile) = @_;
   my $fh;
   if ($infile =~ /\.gz$/){
      open($fh, '-|', "gunzip -c $infile") or die "Could not open $infile $!\n";
   } else {
      open($fh, '<', $infile) or die "Could not open $infile $!\n";
   }
   my %sizes = ();
   while(<$fh>){
      chomp;
      next if /^chrom/;
      my ($chrom, $size) = split(/\t/);
      $sizes{$chrom} = $size;
   }
   close($fh);
   return(%sizes);
}

