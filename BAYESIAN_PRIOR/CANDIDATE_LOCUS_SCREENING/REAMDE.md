# What is does
This workflow takes as input the coordinates of a gene and a list of vcf (and bam) files. It searches for the loci in the gene that have not been called in the vcf.
Next it does a pileup in these loci to search for potential false negative candidates. The user receives a list of potential FN candidates, which can be further manually inspected
(for example in IGV).

# Input
  - referece genome
  - bed file of the gene
  - list of bamfiles (one path per line)
  - list of gzvcf files (one path per line)
* The sample IDs in the bam files and the gzvcf files must be the same and in the same order!

The input must be organized in the followinf format (in the ```config.yaml``` file):
```
reference_genome: /hs37d5.fa
gene_bed: /ALK.bed
snakedir: /DATA_LABELING
bamlist_file: /bamlist_file.txt
gzvcf_list_file: /gzvcf_list_file.txt


resources:
  rule_identify_missing_loci: 100000
  rule_make_bam_readcount: 100000
  rule_process_bam_readcount: 100000
```

# How to run

Before running, the working directory should look like this:
```
.
├── Snakefile.smk
├── bamlist_file.txt
├── config.yaml
├── gzvcf_list_file.txt
├── hs37d5.fa -> /data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/T-CELL-PROJECT-DATA/map_to_37/hs37d5.fa
├── hs37d5.fa.fai -> /data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/T-CELL-PROJECT-DATA/map_to_37/hs37d5.fa.fai
└── scripts
    ├── IdentifyMissedLoci.sh
    ├── make_bam_readcounts.sh
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
├── Snakefile.smk
├── bamlist_file.txt
├── config.yaml
├── gzvcf_list_file.txt
├── hs37d5.fa -> /data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/T-CELL-PROJECT-DATA/map_to_37/hs37d5.fa
├── hs37d5.fa.fai -> /data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/T-CELL-PROJECT-DATA/map_to_37/hs37d5.fa.fai
├── missed_loci
│   ├── ...
│   └── ...
├── pileup_missed_loci_processed
│   ├── ...
│   └── ...
├── pileup_missed_loci_raw
│   ├── ...
│   └── ...
└── scripts
    ├── IdentifyMissedLoci.sh
    ├── make_bam_readcounts.sh
    └── process_bam_readcount.py
```

# Output format

| chrom | pos       | ref  | coverage | ref_count | total_alt_count | A_count | C_count | G_count | N_count |
|-------|-----------|------|----------|-----------|----------------|---------|---------|---------|---------|
| 12    | 25358070  | A    | 1        | 0         | 1              | 0       | 0       | 0       | 0       |
| 12    | 25358070  | A    | 1        | 0         | 1              | 0       | 0       | 0       | 0       |
| 12    | 25358071  | A    | 1        | 0         | 1              | 0       | 0       | 0       | 0       |

