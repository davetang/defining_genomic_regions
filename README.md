Defining genomic regions
========================

Clone this repository:

`git clone https://github.com/davetang/defining_genomic_regions.git`

## To run all scripts

`make`

## Annotate BAM files

```
samtools sort my_file.bam my_file
bedtools2/bin/bedtools bamtobed -i my_file.bam > my_file.bed
cat my_file.bed | wc -l
bedtools2/bin/bedtools intersect -a my_file.bed -b gencode_v19_exon_merged.bed.gz -u | wc -l
bedtools2/bin/bedtools intersect -a my_file.bed -b gencode_v19_intergenic.bed.gz -u | wc -l
bedtools2/bin/bedtools intersect -a my_file.bed -b gencode_v19_intron.bed.gz -u | wc -l
```

## Links

See assoicated blog post: <http://davetang.org/muse/2013/01/18/defining-genomic-regions/>
