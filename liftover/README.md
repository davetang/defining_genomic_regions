## README

Visit the [UCSC Genome Browser Store](https://genome-store.ucsc.edu/products/) and download liftOver after creating an account. It is free for personal and non-profit academic research use.

Download a chain file.

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

Check out the chain file.

```bash
zcat hg19ToHg38.over.chain.gz | head -6
chain 20851231461 chr1 249250621 + 10000 249240621 chr1 248956422 + 10000 248946422 2
167376  50041   80290
40302   253649  288020
1044699 1       2
3716    0       3
1134    4       18
```

The [chain format](https://genome.ucsc.edu/goldenPath/help/chain.html) has an initial header line starts with the keyword `chain`, followed by 11 required attribute values, and ends with a blank line. The attributes include:

* `score` -- chain score
* `tName` -- chromosome (reference sequence)
* `tSize` -- chromosome size (reference sequence)
* `tStrand` -- strand (reference sequence)
* `tStart` -- alignment start position (reference sequence)
* `tEnd` -- alignment end position (reference sequence)
* `qName` -- chromosome (query sequence)
* `qSize` -- chromosome size (query sequence)
* `qStrand` -- strand (query sequence)
* `qStart` -- alignment start position (query sequence)
* `qEnd` -- alignment end position (query sequence)
* `id` -- chain ID

The alignment data lines contain three required attribute values:

* `size` -- the size of the ungapped alignment
* `dt` -- the difference between the end of this block and the beginning of the next block (reference sequence)
* `dq` -- the difference between the end of this block and the beginning of the next block (query sequence)

The block chr1:10000-177376 should liftover to the exact coordinates on hg38.

The `liftOver` tool requires four positional arguments: oldFile map.chain newFile unMapped

```bash
perl -le 'print join("\t", "chr1", 10000, 177376)' > chr1_10000_177376.txt
./liftOver chr1_10000_177376.txt hg19ToHg38.over.chain.gz chr1_10000_177376_hg38.txt chr1_10000_177376_unmapped.txt
cat chr1_10000_177376_hg38.txt
# chr1    10000   177376
```

If we create a BED file with a region that doesn't lift over, the output BED file will be trimmed. (The alignment block ends at chr1:10000:177376, so the 1 bp overhang will be trimmed.)

```bash
perl -le 'print join("\t", "chr1", 10000, 177377)' > chr1_10000_177377.txt
./liftOver chr1_10000_177377.txt hg19ToHg38.over.chain.gz chr1_10000_177377_hg38.txt chr1_10000_177377_unmapped.txt
cat chr1_10000_177377_hg38.txt
# chr1    10000   177376
```

