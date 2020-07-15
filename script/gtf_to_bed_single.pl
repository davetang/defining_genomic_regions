#!/usr/bin/env perl
#
# Output BED file with features of interest for genes with only one transcript model
#

use warnings;
use strict;
use Getopt::Std;

my %opts = ();
getopts('i:f:h:', \%opts);

if ($opts{'h'} ||
    !exists $opts{'i'} ||
    !exists $opts{'f'}
){
   usage();
}

my $my_feature = $opts{'f'};
my %gene_anno = ();

# first read through to tally number of transcripts per gene
my $gtf_file = $opts{'i'};
my $gtf = open_file($gtf_file);
while(<$gtf>){
   chomp;
   next if /^#/;
   my ($sequence, $source, $feature, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/);

   my $gene_id = '.';
   if ($attributes =~ /gene_id\s"([a-zA-Z0-9._]+)";/){
      $gene_id = $1;
   }

   if ($feature eq "transcript"){
      if (exists $gene_anno{$gene_id}->{'COUNT'}){
         $gene_anno{$gene_id}->{'COUNT'}++;
      } else {
         $gene_anno{$gene_id}->{'COUNT'} = 1;
      }
   }

}
close($gtf);

# second read through to output features of interest for genes with only one transcript model
my $gtf2 = open_file($gtf_file);
while(<$gtf2>){
   chomp;
   next if /^#/;
   my ($sequence, $source, $feature, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/);

   my $gene_id = '.';
   if ($attributes =~ /gene_id\s"([a-zA-Z0-9._]+)";/){
      $gene_id = $1;
   }

   if ($gene_anno{$gene_id}->{COUNT} == 1 && $feature eq $my_feature){
      # BED is 0-based
      $start -= 1;
      print join("\t", $sequence, $start, $end, $gene_id, $score, $strand), "\n";
   }

}
close($gtf2);

sub open_file {
   my ($infile) = @_;
   my $fh;
   if ($infile =~ /\.gz$/){
      open($fh, '-|', "gunzip -c $infile") || die "Could not open $infile $!\n";
   } else {
      open($fh, '<', $infile) || die "Could not open $infile $!\n";
   }
   return($fh);
}


sub usage {
print STDERR <<EOF;
Usage: $0 -f FILE -l STRING

Where:   -i         GTF file
         -f         Feature to output, e.g. gene, transcript, CDS, exon
         -h         this helpful usage message

EOF
exit();
}

__END__

