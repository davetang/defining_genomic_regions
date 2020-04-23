## README

Download GTF file

    wget -c -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz

Convert to BED using "gene" entries in the GTF file

    ./gtf_to_bed.pl -i gencode.v33.annotation.gtf.gz -f gene > gencode.v33.gene.bed

