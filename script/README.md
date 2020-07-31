## README

Download GTF file

    wget -c -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz

Convert to BED using "gene" entries in the GTF file

    ./gtf_to_bed.pl -i gencode.v33.annotation.gtf.gz -f gene > gencode.v33.gene.bed

## Scripts

The `merge_by_id.pl` script will merge by the annotation column (column 4). This is useful when you only want to merge features with the same ID; I could not find an elegant solution, so I wrote this script. It requires `bedtools` and will create a lot of temporary files (one per feature in the BED file). It works by running `bedtools merge` on each feature, so it can be very slow with files with a lot of different features. You can specify more threads if you have multiple cores on your system.

