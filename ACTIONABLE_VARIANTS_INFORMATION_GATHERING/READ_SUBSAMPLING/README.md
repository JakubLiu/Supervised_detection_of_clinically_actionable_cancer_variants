## General information
The goal if this small pipeline is to perform read subsampling and then call somatic variants on the various levels of read subsamplings. By doing that we can asses the effect of coverage on the variants calling output.

## I/O
#### Input
- a .bam file or a list of .bam files
- a reference genome
- desired proportion(s) of reads to keep per .bam file

The input is formatted as a **config** file:

```
# paths (to the bam files) can be absolute or relative
bam_files:
  - data/bams/A.bam
  - data/bams/B.bam
  - bams/C.bam

# percentage of reads to keep
proportion_reads_to_keep:
  - 0.2
  - 0.5
  - 0.7

reference_genome: data/genome.fa
```

#### Output
- bamfiles with subsampled reads (each combination of (.bam x subsampling_proporion))
- VCF files based on the subsampled .bam files

## How to run
#### prerequisites
- the correct conda environment based on: ```save_the_environment.yaml```
- the reference genome needs the following supporting files:
```
genome.fa.sa
genome.fa.pac
genome.fa.fai
genome.fa.bwt
genome.fa.ann
genome.fa.amb
genome.dict
genome.fa --> the reference genome itself
```
```
snakemake --snakefile Snakefile.smk --configfile <your config.yaml file> --cores <the number of cores>
```

or as a Surm job (BIH HPC specific):

```
snakemake --snakefile Snakefile.smk --configfile <your config.yaml file> --profile=cubi-v1 --jobs 1
```

#### (optional)
Often there is the error that stemms from **Mutect2**. This software expects the readgroups to be present in the .bam files.
To solve this issue do the following for each .bam file.
##### 1.) add the readgroups:
```
gatk AddOrReplaceReadGroups \
  -I <your input .bam> \
  -O <name of your output .bam> \
  -RGID 1 \     # this is just a dummy name in most cases
  -RGLB lib1 \  # this is just a dummy name in most cases
  -RGPL illumina \   # sequencing platform
  -RGPU unit1 \     # this is just a dummy name in most cases
  -RGSM <name of your samples (e.g. patient_1)>
```

##### 2.) sort the .bam file (the one with the readgroups)
```
samtools sort <input bam> -o <output bam>
```

##### 3.) index the sorted .bam file
```
samtools index <input bam>
```

An example usage setup can be found in ```/READ_SUBSAMPLING/EXAMPLE_RUN/```
