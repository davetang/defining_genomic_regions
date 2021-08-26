#!/usr/bin/env perl
#
# Some GTF files do not contain exon features. This script add exon features based on the CDS feature.
# This script requires the gene_id attribute and will set the transcript_id to the gene_id.
# Thus only use this script for GTF files where there is only one transcript per gene.
# If you have more than one transcript per gene, remove or comment out line 47.
#
# A CDS is a contiguous sequence which begins with, and includes, the start codon but does not include the stop codon.
# An exon is a region of the transcript sequence within a gene which is not removed from the primary RNA transcript by RNA splicing.
#

use warnings;
use strict;
use Getopt::Std;

my %opts = ();
getopts('i:f:h:u:d:', \%opts);

if ($opts{'h'} ||
    !exists $opts{'i'}
){
   usage();
}

my $gtf = $opts{'i'};

if ($gtf =~ /\.gz$/){
   open(IN, '-|', "gunzip -c $gtf") || die "Could not open $gtf: $!\n";
} else {
   open(IN, '<', $gtf) || die "Could not open $gtf: $!\n";
}

while(<IN>){
   chomp;
   next if /^#/;
   my ($sequence, $source, $feature, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/);

   if ($feature eq 'CDS'){
      my $gene_id = '';
      if ($attributes =~ /gene_id\s"([a-zA-Z0-9._]+)";/){
         $gene_id = $1;
      } else {
         die "Could not extract gene_id on line $.: $_\n";
      }
      print "$_\n";
      $attributes = "gene_id \"$gene_id\"; transcript_id \"$gene_id\";";
      print join("\t", $sequence, $source, 'exon', $start, $end + 3, $score, $strand, $phase, $attributes), "\n";
   } else {
      print "$_\n";
   }

}
close(IN);

sub usage {
print STDERR <<EOF;
Usage: $0 -i FILE

Where:   -i         GTF file
         -h         this helpful usage message

EOF
exit();
}

__END__

