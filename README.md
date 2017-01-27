Defining genomic regions
========================

Clone this repository:

~~~~{.bash}
git clone https://github.com/davetang/defining_genomic_regions.git
~~~~

## To run all scripts

~~~~{.bash}
make

# if everything ran successfully
for file in `ls *.gz`; do md5sum $file; done
bd83e28270e595d3bde6bfcb21c9748f  gencode.v19.annotation.gtf.gz
8c97ec4b54eaa176ba1e48bfeb60c08a  gencode_v19_exon_merged.bed.gz
ea03038b873ba2612383a4c0949c835d  gencode_v19_intergenic.bed.gz
4d5ff850e3115077bf50d87bc406a84f  gencode_v19_intron.bed.gz
48dbe15f4498baad1a2327c774a692c8  promoter.bed.gz
9d513cad3aafd5690bf8bbebb24b4df4  transcript_utr.bed.gz
35aed6aac655182c653cdc72060b914d  transcript_utr_number.out.gz
~~~~

## Annotate BAM files

~~~~{.bash}
samtools sort my_file.bam my_file
bedtools2/bin/bedtools bamtobed -i my_file.bam > my_file.bed
cat my_file.bed | wc -l
bedtools2/bin/bedtools intersect -a my_file.bed -b gencode_v19_exon_merged.bed.gz -u | wc -l
bedtools2/bin/bedtools intersect -a my_file.bed -b gencode_v19_intergenic.bed.gz -u | wc -l
bedtools2/bin/bedtools intersect -a my_file.bed -b gencode_v19_intron.bed.gz -u | wc -l
~~~~

## Links

See assoicated blog post: <https://davetang.org/muse/2013/01/18/defining-genomic-regions/>

