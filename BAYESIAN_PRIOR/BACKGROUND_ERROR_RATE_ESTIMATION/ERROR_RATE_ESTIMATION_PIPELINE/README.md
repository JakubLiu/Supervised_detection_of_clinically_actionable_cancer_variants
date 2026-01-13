# What it does
This workflow runs a background error rate estimation for all loci in a given gene in all samples.
For locus _j_ in sample _i_ and nucleotide _k_ the error rate is defined as:


$$
\hat{e}_{ijk} = X_{ijk}/(n_{ij} + \epsilon)
$$

Where, $X_{ijk}$ is the number of reads supporting nucleotide $k$ at locus $j$ in sample $i$.
$n_{ij}$ is the coverage, and $\epsilon$ is a pseudocount of $10^{-10}$ added to the coverage
to avoid division by zero errors in loci with no coverage.

## Caveats
- results are stable only for SNVs, indels or CNVs are not supported
- non-standard bases (like $N$) are not supported (removed in the scripts)


# Expected inputs
1.) reference genome (and all the index files needed by samtools)
2.) list of bamfiles (paths to the bamfiles, one path per line)
    2.1.) the corresponding .bai files must be in the same directory
3.) gene bedfile
```
chrom  start  stop
```
  3.1.) the chromosome identifier must match the one used in the reference genome
4.) gene name
5.) path to the current directory

All the inputs must be placed in the ```config.yaml``` file that has the following format:
```
reference_genome: hs37d5.fa
bamlist_file: /ERROR_RATE_ESTIMATION_PIPELINE/bamlist_file.txt
bamfile_extension: .dedupa.rg.sorted.bam
gene_bedfile: /ERROR_RATE_ESTIMATION_PIPELINE/AKT1.bed
gene_name: AKT1
snakedir: /ERROR_RATE_ESTIMATION_PIPELINE



# resources:
resources_mb:
  rule_pileup: 100000
  rule_process_pileup: 100000
  rule_estimate_error_rate: 100000
```

# How to run
The directory tree should look like this before running:
```
.
├── AKT1.bed
├── Snakefile.smk
├── bamlist_file.txt
├── config.yaml
└── scripts
    ├── error_rate.py
    ├── pileup.sh
    └── process_pileup.py
```

To run type:
```
snakemake --snakefile Snakefile.smk --configfile config.yaml --cores <num cores>
```
or via Slurm on the BIH HPC:
```
snakemake --snakefile Snakefile.smk --configfile config.yaml --profile=cubi-v1 --jobs <num jobs>
```

After succesfull completion the directory tree should look like this:

```
├── AKT1.bed
├── Snakefile.smk
├── bamlist_file.txt
├── config.yaml
├── error_rate
│   └── AKT1.error_rates.txt
├── pileup
│   └── AKT1.pileup.txt
├── procesed_pileup
│   └── AKT1.pileup.processed.txt
└── scripts
    ├── error_rate.py
    ├── pileup.sh
    └── process_pileup.py

```

# Output
The output has the following form:

| sampleID | chrom | start | stop | ref | coverage | n_alt_bases | read_pileup | $\hat{e}_{ijA}$ | $\hat{e}_{ijT}$ | $\hat{e}_{ijG}$ | $\hat{e}_{ijC}$ |
|---------|---------|---------|---------|---------|---------|---------|---------|---------|----------|----------|----------|
|  Steve's lobster sample 1   |  14 |  100 | 2000 | $A$ | 32| 2| AAAAAA...| 0.00 | 0.00 | 0.001 | 0.00 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |

## Caveat
The error rate for the reference nucleotide is always set to $0.0$. For example:
$ref_{j} = A$ --> $\hat{e}_{ijA} = 0.0$
