## README

Setup latest version of BEDTools.

```bash
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz
tar -xzf bedtools-2.29.2.tar.gz
cd bedtools2
make all

cd ~/bin/
ln -s ~/src/bedtools2/bin/bedtools

bedtools --version
# bedtools v2.29.2
```

Download GENCODE GTF files.

```bash
parallel --verbose wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{}/gencode.v{}.annotation.gtf.gz ::: {20..35}
```

Calculate stats.

```bash
parallel --verbose "../run.pl {} ../chrom_info/hg38.genome > {.}.stats" ::: *.gtf.gz
```

Plot using `plot_stats.Rmd`.

<img src="https://github.com/davetang/defining_genomic_regions/blob/main/gencode/genomic_region_proportion.png" width="600" />
<img src="https://github.com/davetang/defining_genomic_regions/blob/main/gencode/genomic_region_length.png" width="600" />
<img src="https://github.com/davetang/defining_genomic_regions/blob/main/gencode/gene_type.png" width="600" />

