#!/bin/env perl

use strict;
use warnings;

my $usage = "Usage: $0 <infile.annotation.gtf> <padding>\n";
my $infile = shift or die $usage;
my $span = shift or die $usage;

if ($span !~ /^\d+$/){
   die "Please enter a numeric value for the padding\n";
}

my $hg19 = 'hg19.genome';
my %hg19 = ();

open(IN,'<',$hg19) || die "Could not open $hg19: $!\n";
while(<IN>){
   chomp;
   #chr9_gl000201_random    36148
   my ($chr, $end) = split(/\t/);
   $hg19{$chr} = $end;
}
close(IN);

if ($infile =~ /\.gz/){
   open(IN,'-|',"gunzip -c $infile") || die "Could not open $infile: $!\n";
} else {
   open(IN,'<',$infile) || die "Could not open $infile: $!\n";
}

while(<IN>){
   chomp;
   next if (/^#/);
   #chr11   HAVANA  transcript      65265233        65273940        .       +       .       gene_id "ENSG00000251562.3"; transcript_id "ENST00000534336.1"; gene_type "processed_transcript"; gene_status "KNOWN"; gene_name "MALAT1"; transcript_type "non_coding"; transcript_status "KNOWN"; transcript_name "MALAT1-001"; level 2; havana_gene "OTTHUMG00000166322.1"; havana_transcript "OTTHUMT00000389143.1";
   my ($chr,$source,$type,$start,$end,$score,$strand,$phase,$annotation) = split(/\t/);
   next unless $type eq 'transcript';
   my @annotation = split(/;\s/,$annotation);
   my $transcript_id = 'none';
   foreach my $blah (@annotation){
      my ($type,$name) = split(/\s+/,$blah);
      if ($type eq 'transcript_id'){
         $transcript_id = $name;
         $transcript_id =~ s/"//g;
      }
   }
   if ($transcript_id eq 'none'){
      die "No name for entry $.\n";
   }
   my $promoter_start = '';
   my $promoter_end = '';
   if ($strand eq '+'){
      $promoter_start = $start - $span;
      $promoter_end = $start + $span;
   } else {
      $promoter_start = $end - $span;
      $promoter_end = $end + $span;
   }
   if ($promoter_start < 0){
      warn "Adjusted promoter start to 0\n";
      $promoter_start = 0;
   } elsif ($promoter_end > $hg19{$chr}){
      warn "Adjusted promoter end to $hg19{$chr}\n";
      $promoter_end = $hg19{$chr};
   }
   print join("\t",$chr,$promoter_start,$promoter_end,$transcript_id,0,$strand),"\n";
}
close(IN);

exit(0);
