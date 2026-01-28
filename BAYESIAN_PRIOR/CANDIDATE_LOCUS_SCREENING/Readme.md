# What is does
This workflow takes as input the coordinates of a gene and a list of vcf (and bam) files. It searches for the loci in the gene that have not been called in the vcf.
Next it does a pileup in these loci to search for potential false negative candidates. The user receives a list of potential FN candidates, which can be further manually inspected
(for example in IGV).

# Input
  - referece genome
  - bed file of the gene
  - list of bamfiles (one path per line) one for the tumor and one for the normal sample
  - list of gzvcf files (one path per line)
* The sample IDs in the bam files and the gzvcf files must be the same and in the same order!

The input must be organized in the followinf format (in the ```config.yaml``` file):
```
reference_genome: /DATA_LABELING/hs37d5.fa
gene_bed: /DATA_LABEL_POLIGON/ALK.bed
snakedir: /DATA_LABELING
bamlist_file_tumor: /DATA_LABELING/bamlist_file_tumor.txt
bamlist_file_normal: /DATA_LABELING/bamlist_file_normal.txt
gzvcf_list_file: /DATA_LABELING/gzvcf_list_file.txt


resources:
  rule_identify_missing_loci: 100000
  rule_make_bam_readcount_tumor: 100000
  rule_process_bam_readcount_tumor: 100000
  rule_make_bam_readcount_normal: 100000
  rule_process_bam_readcount_normal: 100000
  rule_merge: 100000



```

# How to run

Before running, the working directory should look like this:
```
.
├── Snakefile.smk
├── Untitled.ipynb
├── bamlist_file_normal.txt
├── bamlist_file_tumor.txt
├── config.yaml
├── gzvcf_list_file.txt
├── hs37d5.fa
├── hs37d5.fa.fai
└── scripts
    ├── IdentifyMissedLoci.sh
    ├── make_bam_readcounts.sh
    ├── merge.py
    └── process_bam_readcount.py

```

Run like this:

```
snakemake --snakefile Snakefile.smk --configfile <the config.yaml file> --cores <number of cores>
```

Or on the BIH HPC

```
snakemake --snakefile Snakefile.smk --configfile <the config.yaml file> --profile=cubi-v1 --jobs <num jobs>
```

After running the working directory should look like this:

```
.
├── Snakefile.smk
├── Untitled.ipynb
├── bamlist_file_normal.txt
├── bamlist_file_tumor.txt
├── config.yaml
├── gzvcf_list_file.txt
├── hs37d5.fa
├── hs37d5.fa.fai
├── merged
│   ├── ...
│   └── ...
├── missed_loci
│   ├── ...
│   └── ...
├── pileup_missed_loci_processed
│   ├── ...
│   └── ...
├── pileup_missed_loci_processed_normal
│   ├── ...
│   └── ...
├── pileup_missed_loci_raw
│   ├── ...
│   └── ...
├── pileup_missed_loci_raw_normal
│   ├── ...
│   └── ...
└── scripts
    ├── IdentifyMissedLoci.sh
    ├── make_bam_readcounts.sh
    ├── merge.py
    └── process_bam_readcount.py
```

# Output format

| chrom | pos | ref | coverage_x | ref_count_x | total_alt_count_x | A_count_x | C_count_x | G_count_x | N_count_x | coverage_y | ref_count_y | total_alt_count_y | A_count_y | C_count_y | G_count_y | N_count_y |
|-------|-----|-----|------------|-------------|-------------------|-----------|-----------|-----------|-----------|------------|-------------|-------------------|-----------|-----------|-----------|-----------|
| 2 | 29415957 | T | 21 | 20 | 1 | 0 | 1 | 0 | 0 | 7 | 7 | 0 | 0 | 0 | 0 | 0 |
| 2 | 29415957 | T | 21 | 20 | 1 | 0 | 1 | 0 | 0 | 7 | 7 | 0 | 0 | 0 | 0 | 0 |
| 2 | 29415957 | T | 21 | 20 | 1 | 0 | 1 | 0 | 0 | 7 | 7 | 0 | 0 | 0 | 0 | 0 |
| 2 | 29415957 | T | 21 | 20 | 1 | 0 | 1 | 0 | 0 | 7 | 7 | 0 | 0 | 0 | 0 | 0 |



