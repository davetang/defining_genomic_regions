## README

Download genome annotations.

```bash
for genome in hg19 hg38; do
   mysql --user=genome \
         --host=genome-mysql.cse.ucsc.edu \
         -A \
         -e "select chrom, size from ${genome}.chromInfo" | grep -v "^chrom" > ${genome}.genome
done

wget -q -O - http://genome-test.cse.ucsc.edu/~hiram/hubs/Plants/araTha1/araTha1.chrom.sizes |
  sed 's/^chr//' |
  sed 's/Cp/Pt/' > araTha1.genome
```

