## General outline
This pipeline aims to quantify the amount of clinically actionable variants that have been missed by unbiased and genome wide variants callers (here Mutect2), as
a function of the read coverage depth.

## Description of the rules
#### subsample_bams
- **input**: .bam files (sorted, indexed and with read groups)
- **output**: these .bam files with a reduced read coverage depth
- **what is does**: It uses samtools to keep only a given subset of the reads from the input .bam files

#### mutect2_call & mutect2_filter
- **input**: subsample .bam files (from rule **subsample_bams**)
- **output**: VCF files
- **what they do**: These two rules perform variant calling and filtration of the initially called variants.

#### keep_PASS_variant_calls
- **input**: filtered VCF files
- **output**: VCF files that only include the variant calls that have passed the filters
- **what is does**: It excludes all variants that have failed some of the filters applied by ```gatk FilterMutectCalls```

#### compile_multiple_vcfs
- **input**: a list of VCF files
- **output**: a .txt file showing in how many of these VCF files each variant has been called
- **what is does**: For each unique variant across all VCFs show in which of these files it has been called in.

#### n_out_of_m_samples
- **input**: the output of rule **compile_multiple_vcfs**
- **output**: a .bed file of the coordinates of the variants that have been **NOT CALLED** in n out of m VCFs, where m is the total number of VCFs/patients/samples
- **what is does**: Let's say we have 5 total VCFs (one VCF per patient) and we want to find the coordinates of the variants that have been called in 4/5 samples which is equivalent to bein missed in 1 sample. This rule finds these coordinates.

#### merge_with_actionables_bed
- **input**: a .bed file from rule **n_out_of_m_samples** and an external .bed file of variants of interest (for example clinically actionable variants)
- **output**: a combined .bed file
- **what is does**: A simple horizontal stacking of two files is performed. No check for overlapping or identical regions is performed.

#### mpileup
- **input**: the merged .bed files from rule **merge_with_actionables_bed**, a .bam file and a reference genome (with all it's supporting data structures required by samtools)
- **output**: a samtools mpileup object
- **what is does**: It just runs the samtools mpileup command.

#### postprocess_mpileup
- **input**: the mpileup object from rule **mpileup**
- **output**: a table containing the processed mpileup output and an estimation of whether a given candidate is a False Negative or not (for an example see **EXAMPLE_RUN/postprocessed_mpileup/**)
- **what is does**: It translates the putput of samtools mpileup into actual base sequences. Then based on the locus-specific coverage and the given cutoff value if calls a candidate as a False Negative or not. [1]

[1] Let's say that at our candidate locus we have a read coverage of 100 and we have a threshold of 0.1. If and only if we have at least 0.1*100=10 reads supporting a non-reference allele (at that locus), then the candidate is called a False Negative.



## How to run
#### required inputs
- .bam files (indexed, sorted and with read groups (see gatk AddOrReplaceReadGroups))
- a reference genome along with all the supporting data structures required by Samtools and Mutect2
- the ```SOURCE_DIR/``` that holds the following files:
```
SOURCE_DIR/
  |___ Snakefile.smk
  |___ scripts/
          |__ "all the scripts called by the snakefile"
```
- the conda environment (see **save_the_environment.yaml**)

The inputs should be stated as below in the config file.

```
source_dir: SOURCE_DIR/

bamfile_list:
  - data/samples/A.bam
  - data/samples/B.bam
  - data/samples/C.bam

reference_genome: data/genome.fa

actionable_variants_bed: data/actionable_variants.bed

cutoffs:
  - 0.01    # For example if 0.01 of the reads at a given locus support a non-reference allele, call that candidate as a False Negative
  - 0.05

coverages:
  - 0.25    # keep 0.25 of the reads
  - 0.50

sample_proportions:  # for this we need to know how many samples we have in total
  - 1        # not called in 1 out of the total number of samples
  - 2        # not called in 2 out of the total number of samples
```

If you wish to keep all the intermediate directories and files, run:

```
snakemake --snakefile SOURCE_DIR/Snakefile.smk --configfile <your config yaml file> --config remove_intermediate_dirs=False --cores <number of cores>
```

If you with all files and directories to be removed once they become unnecessary, run:

```
snakemake --snakefile SOURCE_DIR/Snakefile.smk --configfile <your config yaml file> --cores <number of cores>
```

How to run on the BIH HPC cluster:
  - run it on a detachable linux screen and prefferably on a login node
  - especially if you run the pipeline on an HPC as a Slurm job, it is important to reference the SOURCE_DIR/, else Snakemake might not find the scripts/ directory
```
snakemake --snakefile SOURCE_DIR/Snakefile.smk --configfile <your config yaml file> --profile=cubi-v1 --jobs 1
```
