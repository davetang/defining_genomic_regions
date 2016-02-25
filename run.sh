#!/bin/bash

echo Checking for bedtools

if [ ! -d bedtools2 ]
then
   git clone https://github.com/arq5x/bedtools2.git
   cd bedtools2
   make clean; make all
   cd ..
fi

echo Downloading GENCODE annotations

v=19

if [ ! -f gencode.v$v.annotation.gtf.gz ]
then
   wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_$v/gencode.v$v.annotation.gtf.gz
fi

echo Creating exonic regions

if [ ! -f gencode_v${v}_exon_merged.bed.gz ]
then
   gunzip -c gencode.v$v.annotation.gtf.gz |
   awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' |
   bedtools2/bin/sortBed |
   bedtools2/bin/mergeBed -i - | gzip > gencode_v${v}_exon_merged.bed.gz
fi

echo Creating intronic regions

if [ ! -f gencode_v${v}_intron.bed.gz ]
then
   gunzip -c gencode.v$v.annotation.gtf.gz |
   awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |
   bedtools2/bin/sortBed |
   bedtools2/bin/subtractBed -a stdin -b gencode_v${v}_exon_merged.bed.gz |
   gzip > gencode_v${v}_intron.bed.gz
fi

# echo Downloading hg19 coordinates
# 
# if [ ! -f hg19.genome ]
# then
#    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
#    "select chrom, size from hg19.chromInfo"  > hg19.genome
# fi

echo Creating intergenic regions

if [ ! -f gencode_v${v}_intergenic.bed.gz ]
then 
   gunzip -c gencode.v$v.annotation.gtf.gz |
   awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |
   sort -k1V -k2,2n |
   bedtools2/bin/complementBed -i stdin -g hg19.genome |
   gzip > gencode_v${v}_intergenic.bed.gz
fi

echo Counting UTRs

if [ ! -f transcript_utr_number.out.gz ]
then
   perl check_utr.pl gencode.v19.annotation.gtf.gz | gzip > transcript_utr_number.out.gz
fi

echo Creating UTRs

if [ ! -f transcript_utr.bed.gz ]
then
   perl print_utr.pl gencode.v19.annotation.gtf.gz | gzip > transcript_utr.bed.gz
fi

echo Creating promoter region

if [ ! -f promoter.bed.gz ]
then
   perl promoter.pl gencode.v19.annotation.gtf.gz 200 | gzip > promoter.bed.gz
fi


echo Done
