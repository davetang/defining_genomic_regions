## README

Guess the genome assembly from an unknown BED file by checking whether coordinates span outside the defined boundary and checking the overlap with the expected genomics feature. For example, I have two BED files containing coordinates for exons but I am not sure what genome assembly the coordinates coorespond to.

## Boundary check

Human genome chromosome sizes can be downloaded from the UCSC Genome Browser's database. I have included the two genome size files in this repository, so you do not need to run the command below. If you want to run it yourself, you will need Docker.

```bash
docker run --rm -u $(stat -c "%u:%g" $HOME) -v $(pwd):$(pwd) -w $(pwd) mariadb:10.3 mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -P 3306 -e "select chrom, size from hg19.chromInfo" | gzip > hg19.genome.gz
docker run --rm -u $(stat -c "%u:%g" $HOME) -v $(pwd):$(pwd) -w $(pwd) mariadb:10.3 mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -P 3306 -e "select chrom, size from hg38.chromInfo" | gzip > hg38.genome.gz
```

Use `check_size.pl` to check.

```bash
./check_size.pl -s hg19.genome.gz -b unknown1.bed.gz # lots of warnings
./check_size.pl -s hg38.genome.gz -b unknown1.bed.gz # no warnings

./check_size.pl -s hg19.genome.gz -b unknown2.bed.gz # no warnings
./check_size.pl -s hg38.genome.gz -b unknown2.bed.gz # lots of warnings
```

## Overlap check

We will use the UCSC Genome Browser's database and the RefSeq database to create a BED file containing exonic regions. UCSC Genome Browser's [internal database representations](https://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1) of coordinates always have a zero-based start and a one-based end, so we do not need to change the coordinates.

```bash
docker run --rm -u $(stat -c "%u:%g" $HOME) -v $(pwd):$(pwd) -w $(pwd) mariadb:10.3 mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -P 3306 -D hg19 -e 'select chrom,exonStarts,exonEnds from refGene' | split_exon.pl | grep -v "_" | sort -k1,1V -k2,2n | uniq | gzip > hg19.genes.gz
docker run --rm -u $(stat -c "%u:%g" $HOME) -v $(pwd):$(pwd) -w $(pwd) mariadb:10.3 mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -P 3306 -D hg38 -e 'select chrom,exonStarts,exonEnds from refGene' | split_exon.pl | grep -v "_" | sort -k1,1V -k2,2n | uniq | gzip > hg38.genes.gz
```

Use `bedtools jaccard` to calculate the [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index).

```bash
# low Jaccard index
bedtools jaccard -a unknown1.bed.gz -b hg19.genes.gz | column -t
intersection  union      jaccard    n_intersections
3942923       135026099  0.0292012  25031

# much higher Jaccard index
bedtools jaccard -a unknown1.bed.gz -b hg38.genes.gz | column -t
intersection  union      jaccard  n_intersections
36268855      103033712  0.35201  190064

bedtools jaccard -a unknown2.bed.gz -b hg38.genes.gz | column -t
intersection  union      jaccard    n_intersections
3566997       130829021  0.0272646  24223

bedtools jaccard -a unknown2.bed.gz -b hg19.genes.gz | column -t
intersection  union      jaccard  n_intersections
32064945      101997528  0.31437  186043
```

Coordinates for `unknown1.bed.gz` are probably for hg38 and `unknown2.bed.gz` are for hg19.

## Padding check

In addition I want to check whether the coordinates are "padded", which means that additional bps are added. If coordinates are padded, then if I shorten regions, the Jaccard index should increase.

```bash
# first check size of smallest region
zcat unknown1.bed.gz | perl -lane 'print $F[2] - $F[1]' | sort -n | head -1
110

# remove 50 bps from start and end
zcat unknown1.bed.gz | perl -lane 'print join("\t", $F[0], $F[1]+50, $F[2]-50)' | gzip > unknown1_shortened.bed.gz

bedtools jaccard -a unknown1_shortened.bed.gz -b hg38.genes.gz | column -t
intersection  union     jaccard   n_intersections
27250690      92010876  0.296168  189508
```

The Jaccard index is decreased, so the BED file is probably not padded.

If the original BED file was padded say by 50 bp, we would have a Jaccard index of 0.31 and "removing" the padding would increase the Jaccard index to 0.35.

```bash
zcat unknown1.bed.gz | perl -lane 'print join("\t", $F[0], $F[1]-50, $F[2]+50)' | gzip > unknown1_lengthened.bed.gz
bedtools jaccard -a unknown1_lengthened.bed.gz -b hg38.genes.gz | column -t
intersection  union      jaccard   n_intersections
37895511      120602616  0.314218  190090
```

