#!/usr/bin/env perl
#
# Some GTF files do not contain exon features and thus do not work well with some tools. This script add exon features based on the CDS feature.
#
# Furthermore, some GTF files set all transcript_id's to "unknown_transcript_1", which can also create problems.
#
# This script can replace transcript_id's with the gene_id. However, this creates a problem when there are more than one transcript per gene.
# Thus only use the -t option when there is only one transcript model per gene.
#
# Lastly, the CDS, start_codon, and stop_codon lines will not be outputted by default. Use the -c option to output them.
#
# For your information, below are the definitions of CDS and exon:
#
# A CDS is a contiguous sequence which begins with, and includes, the start codon but does not include the stop codon.
# An exon is a region of the transcript sequence within a gene which is not removed from the primary RNA transcript by RNA splicing.
#

use warnings;
use strict;
use Getopt::Std;

my %opts = ();
getopts('i:h:c:t:', \%opts);

if ($opts{'h'} ||
    !exists $opts{'i'}
){
   usage();
}

my $gtf = $opts{'i'};
my $keep_cds = 0;
my $replace_tid = 0;

if (exists $opts{'c'}){
   $keep_cds = 1;
}
if (exists $opts{'t'}){
   $replace_tid = 1;
}

# store coordinates to check for overlaps
my %cds = ();
my %exon = ();

if ($gtf =~ /\.gz$/){
   open(IN, '-|', "gunzip -c $gtf") || die "Could not open $gtf: $!\n";
} else {
   open(IN, '<', $gtf) || die "Could not open $gtf: $!\n";
}

LINE: while(<IN>){
   chomp;
   next if /^#/;
   my ($sequence, $source, $feature, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/);

   my $gene_id = '';
   if ($attributes =~ /gene_id\s"([\/a-zA-Z0-9._-]*)";/){
      $gene_id = $1;
   } else {
      die "[ERROR] Could not extract gene_id on line $.: $_\n";
   }
   if ($gene_id eq ''){
      warn("[WARNING] $feature on line $. is not associated with any gene_id: skipping\n");
      warn("[WARNING] $_\n");
      next LINE;
   }

   if ($replace_tid){
      $attributes =~ s/transcript_id "[\/a-zA-Z0-9._-]+"/transcript_id "$gene_id"/;
   }

   if ($feature eq 'CDS'){
      if ($keep_cds){
         print "$_\n";
      }
      print join("\t", $sequence, $source, 'exon', $start, $end + 3, $score, $strand, $phase, $attributes), "\n";
      $cds{$start} = $end + 3;
   } elsif ($feature eq 'start_codon') {
      if ($keep_cds){
         print join("\t", $sequence, $source, $feature, $start, $end, $score, $strand, $phase, $attributes), "\n";
      }
   } elsif ($feature eq 'stop_codon') {
      if ($keep_cds){
         print join("\t", $sequence, $source, $feature, $start, $end, $score, $strand, $phase, $attributes), "\n";
      }
   } elsif ($feature eq 'exon') {
      $exon{$start} = $end;
   } else {
      print join("\t", $sequence, $source, $feature, $start, $end, $score, $strand, $phase, $attributes), "\n";
   }

}
close(IN);

# Issue warning if the newly created exon has coordinates identical to an existing exon
foreach my $start (keys %cds){
   my $end = $cds{$start};
   if (exists $exon{$start} && $exon{$start} == $end){
      warn("[WARNING] Exon $start-$end is repeated; please confirm that they below to different transcript_id's.\n");
   }
}

warn("[WARNING] Finished processing $gtf.\n");

sub usage {
print STDERR <<EOF;
Usage: $0 -i FILE

Where:   -i         GTF file
         -c         output CDS, start_codon, and stop_codon (default: FALSE)
         -t         replace transcript ID with gene ID (default: FALSE)
         -h         this helpful usage message

EOF
exit();
}

__END__

