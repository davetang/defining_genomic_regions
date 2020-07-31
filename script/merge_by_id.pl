#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;
use File::Path;

my %opts = ();
getopts('h:f:t:', \%opts);

if ($opts{'h'} ||
    !exists $opts{'f'}
){
   usage();
}

chomp(my $bedtools = `command -v bedtools`);
if ($bedtools eq ''){
   die "Could not find bedtools\n";
}

my $infile = $opts{'f'};
my $fork_process = 1;
if (exists $opts{'t'}){
   $fork_process = $opts{'t'};
}
warn("Using $fork_process threads\n");

my $fh;
if ($infile =~ /\.gz$/){
   open($fh, '-|', "gunzip -c $infile") || die "Could not open $infile: $!\n";
} else {
   open($fh, '<', $infile) || die "Could not open $infile: $!\n";
}

# store all IDs
my %all_id = ();
while(<$fh>){
   chomp;
   my ($chr, $start, $end, $id, @rest) = split(/\t/);
   $all_id{$id} = 1;
}
close($fh);

my $tmp_dir = time() . "_tmp";
mkdir($tmp_dir) || die "Could not create $tmp_dir: $!\n";

my @command = ();
foreach my $id (keys %all_id){
   my $command;
   if ($infile =~ /\.gz$/){
      $command = "gunzip -c $infile | grep $id | sort -k1,1V -k2,2n | $bedtools merge -i - > $tmp_dir/$id.bed";
      push(@command, $command);
   } else {
      $command = "cat $infile | grep $id | sort -k1,1V -k2,2n | $bedtools merge -i - > $tmp_dir/$id.bed";
      push(@command, $command);
   }
}

my @child = ();
while(scalar(@command) > 0){
   for (1 .. $fork_process){
      my $pid = fork();
      if ($pid) {
         # parent
         push(@child, $pid);
         pop(@command);
      } elsif ($pid == 0) {
         # child
         if (scalar(@command) > 0){
            # print "$command[-1]\n";
            system($command[-1]);
         }
         exit(0);
      } else {
         die "Couldn't fork: $!\n";
      }
   }
   foreach my $pid (@child) {
      waitpid($pid, 0);
   }
}

# oepn merged filed
opendir(DIR, $tmp_dir) || die "Could not open $tmp_dir: $!\n";
while(my $bed = readdir(DIR)){
   next unless $bed =~ /\.bed$/;
   my $id = $bed;
   $id =~ s/\.bed//;
   open(my $fh, '<', "$tmp_dir/$bed") || die "Could not open $tmp_dir/$bed: $!\n";
   while(<$fh>){
      chomp;
      print "$_\t$id\n";
   }
   close($fh);
}
closedir(DIR);

rmtree($tmp_dir) || die "Could not remove $tmp_dir: $!\n";
warn("Done\n");
exit(0);

sub usage {
print STDERR <<EOF;
Usage: $0 -f file -t 16

Where:   -f         BED file with IDs in the fourth column
         -t         threads to use (default 1)
         -h         this helpful usage message

EOF
exit(1);
}

