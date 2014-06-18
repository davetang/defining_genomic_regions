#!/bin/env perl

use strict;
use warnings;

my $usage = "Usage: $0 <infile.annotation.gtf.gz>\n";
my $infile = shift or die $usage;

if ($infile =~ /\.gz/){
   open(IN,'-|',"gunzip -c $infile") || die "Could not open $infile: $!\n";
} else {
   open(IN,'<',$infile) || die "Could not open $infile: $!\n";
}

my %transcript = ();
my $current_transcript = '';

while(<IN>){
   chomp;
   next if (/^#/);
   #chr11   HAVANA  transcript      65265233        65273940        .       +       .       gene_id "ENSG00000251562.3"; transcript_id "ENST00000534336.1"; gene_type "processed_transcript"; gene_status "KNOWN"; gene_name "MALAT1"; transcript_type "non_coding"; transcript_status "KNOWN"; transcript_name "MALAT1-001"; level 2; havana_gene "OTTHUMG00000166322.1"; havana_transcript "OTTHUMT00000389143.1";
   my ($chr,$source,$type,$start,$end,$score,$strand,$phase,$annotation) = split(/\t/);

   my @annotation = split(/;\s/,$annotation);
   my $transcript_id = 'none';

   if ($type eq 'transcript'){
      foreach my $blah (@annotation){
         my ($type,$name) = split(/\s+/,$blah);
         if ($type eq 'transcript_id'){
            $current_transcript = $name;
            $current_transcript =~ s/"//g;
            $transcript{$current_transcript} = 0;
         }
      }
      if ($current_transcript eq 'none'){
         die "No name for entry $.\n";
      }
   }

   if ($type eq 'UTR'){
      $transcript{$current_transcript}++;
   }
}
close(IN);

foreach my $transcript (keys %transcript){
   print "$transcript\t$transcript{$transcript}\n";
}

exit(0);
