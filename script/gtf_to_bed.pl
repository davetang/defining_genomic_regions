#!/usr/bin/env perl
#
# Converts a GTF file into a BED file; the script will use "gene_id" as the BED name, if it exists in the attributes
#

use warnings;
use strict;
use Getopt::Std;

my %opts = ();
getopts('i:f:h:', \%opts);

if ($opts{'h'} ||
    !exists $opts{'f'} ||
    !exists $opts{'i'}
){
   usage();
}

my $gtf = $opts{'i'};
my $my_feature = $opts{'f'};

if ($gtf =~ /\.gz$/){
   open(IN, '-|', "gunzip -c $gtf") || die "Could not open $gtf: $!\n";
} else {
   open(IN, '<', $gtf) || die "Could not open $gtf: $!\n";
}

while(<IN>){
   chomp;
   next if /^#/;
   my ($sequence, $source, $feature, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/);
   next unless $feature eq $my_feature;

   my $name = '.';
   if ($attributes =~ /gene_id\s"([a-zA-Z0-9.]+)";/){
      $name = $1;
   }

   # BED is 0-based
   $start -= 1;
   print join("\t", $sequence, $start, $end, $name, $score, $strand), "\n";
}
close(IN);

sub usage {
print STDERR <<EOF;
Usage: $0 -f FILE -l STRING

Where:   -i         GTF file
         -f         Feature to keep, e.g. gene, transcript, CDS, exon
         -h         this helpful usage message

EOF
exit();
}

__END__

