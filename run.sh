#!/bin/bash

if [ ! -d bedtools2 ]
then
   git clone https://github.com/arq5x/bedtools2.git
   cd bedtools2
   make clean; make all
   cd ..
fi

v=19

if [ ! -f gencode.v$v.annotation.gtf.gz ]
then
   wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_$v/gencode.v$v.annotation.gtf.gz
fi

if [ ! -f gencode_v${v}_exon_merged.bed.gz ]
then
   zcat gencode.v$v.annotation.gtf.gz |
   awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' |
   bedtools2/bin/sortBed |
   bedtools2/bin/mergeBed -i - | gzip > gencode_v${v}_exon_merged.bed.gz
fi

if [ ! -f gencode_v${v}_intron.bed.gz ]
then
   zcat gencode.v$v.annotation.gtf.gz |
   awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |
   bedtools2/bin/sortBed |
   bedtools2/bin/subtractBed -a stdin -b gencode_v${v}_exon_merged.bed.gz |
   gzip > gencode_v${v}_intron.bed.gz
fi

if [ ! -f hg19.genome ]
then
   mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
   "select chrom, size from hg19.chromInfo"  > hg19.genome
fi

if [ ! -f gencode_v${v}_intergenic.bed.gz ]
then 
   zcat gencode.v$v.annotation.gtf.gz |
   awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |
   bedtools2/bin/sortBed |
   bedtools2/bin/complementBed -i stdin -g hg19.genome |
   gzip > gencode_v${v}_intergenic.bed.gz
fi

if [ ! -f transcript_utr_number.out.gz ]
then
   perl check_utr.pl gencode.v19.annotation.gtf.gz | gzip > transcript_utr_number.out.gz
fi

if [ ! -f transcript_utr.bed.gz ]
then
   perl print_utr.pl gencode.v19.annotation.gtf.gz | gzip > transcript_utr.bed.gz
fi

echo Done
