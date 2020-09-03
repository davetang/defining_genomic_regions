#!/usr/bin/env perl
#
# Use script to create exonic, intronic, and intergenic regions from a GTF file
# You must provide a "chrom.size" file, which is simply a flat file that specifies the size of chromosomes
# For an example, take a look at hg19.genome, which was created using the UCSC Genome Browser's MySQL database by running:
#
#    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
#    "select chrom, size from hg19.chromInfo" > hg19.genome
# 

use strict;
use warnings;
use File::Which;
use File::Basename;

my $bedtools = which('bedtools');

if (!$bedtools){
   print STDERR "bedtools was not found in your path:\n\n$ENV{PATH}\n\nPlease install bedtools and add it to your path:\n\n";
   print STDERR "git clone https://github.com/arq5x/bedtools2.git\ncd bedtools2\nmake clean\nmake all\n\n";
   exit(1);
}

my $usage = "Usage: $0 <infile.gtf.gz> <chrom.size>\n";
my $infile = shift or die $usage;
my $genome = shift or die $usage;

if ($infile !~ /\.gtf\.gz$/){
   print STDERR "Please provide a gzipped GTF file\n";
   exit(1);
}

my $basename = basename($infile,  ".gtf.gz");
my $exon_file       = "$basename.exon.merged.bed.gz";
my $intron_file     = "$basename.intron.bed.gz";
my $intergenic_file = "$basename.intergenic.bed.gz";

if (!-e $exon_file){
   warn "Creating exonic regions\n";
   my $command = "gunzip -c $infile | awk 'BEGIN{OFS=\"\\t\";} \$3==\"exon\" {print \$1,\$4-1,\$5}' | bedtools sort | bedtools merge -i - | gzip > $exon_file";
   system($command);
} else {
   warn "$exon_file already exists; skipping exon step\n";
}

if (!-e $intron_file){
   warn "Creating intronic regions\n";
   my $command = "gunzip -c $infile | awk 'BEGIN{OFS=\"\\t\";} \$3==\"gene\" {print \$1,\$4-1,\$5}' | bedtools sort | bedtools subtract -a stdin -b $basename.exon.merged.bed.gz | gzip > $intron_file";
   system($command);
} else {
   warn "$intron_file already exists; skipping intron step\n";
}

if (!-e $intergenic_file){
   warn "Creating intergenic regions\n";
   my $command = "gunzip -c $infile | awk 'BEGIN{OFS=\"\\t\";} \$3==\"gene\" {print \$1,\$4-1,\$5}' | bedtools sort -g $genome | bedtools complement -i stdin -g $genome | gzip > $intergenic_file";
   system($command);
} else {
   warn "$intergenic_file already exists; skipping intergenic step\n";
}

if (-e $exon_file && -e $intron_file && -e $intergenic_file){
   my ($exon_average, $exon_coverage) = stats($exon_file);
   my ($intergenic_average, $intergenic_coverage) = stats($intergenic_file);
   my ($intron_average, $intron_coverage) = stats($intron_file);

   my $total = $exon_coverage + $intergenic_coverage + $intron_coverage;

   printf "exon_coverage: %.2f\n", $exon_coverage*100/$total;
   printf "intron_coverage: %.2f\n", $intron_coverage*100/$total;
   printf "intergenic_coverage: %.2f\n", $intergenic_coverage*100/$total;
   print "exon_length: $exon_average\n";
   print "intron_length: $intron_average\n";
   print "intergenic_length: $intergenic_average\n";
}

sub stats {

   my ($infile) = @_;
   my $coverage = 0;
   my $total = 0;
   my $average = 0;

   open(IN, '-|' ,"gunzip -c $infile") || die "Could not open $infile: $!\n";
   while(<IN>){
      chomp;
      ++$total;
      my ($chr, $start, $end) = split(/\t/);
      my $c = $end - $start;
      $coverage += $c;
   }
   close(IN);

   $average = sprintf("%.2f", $coverage / $total);

   return($average, $coverage);
}

exit(0);

