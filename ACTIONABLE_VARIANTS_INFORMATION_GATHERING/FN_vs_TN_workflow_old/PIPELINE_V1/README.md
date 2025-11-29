## Theoretical background

The goal of this PhD project is to increase the senistivity of variant calling of actionable variants. Therefore loci that are not being called by existing software are of interest.
If a given locus is not being called, it's due to these two main reasons:
  - it is truly biologically absent (no variant) --> true negative (TN)
  - it is a true variant that is missed by the software --> false negative (FN) --> **we want to capture these**

To achieve this a mathematical model is needed (most probably a Deep Learning model). Therefore training data, where the targets will be the TN and FN sites, is needed.
Since there exists no ground truth for this data, we need to guess whether a given locus that is not called is a TN or an FN. This is the task of the following workflow.

We have data from patients with multiple tumor biopsies.

<p align="center">
  <img src="https://github.com/JakubLiu/Supervised_detection_of_clinically_actionable_cancer_variants/blob/main/ACTIONABLE_VARIANTS_INFORMATION_GATHERING/FN_vs_TN_workflow/venn2.png.png" width="300">
</p>

If a variant has been not consistently called across all biopsy samples from the same patient (e.g. called in 2/3 biopsies), then it is of interest to us
and we will reffer to it as a **candidate**.
We would like to know if it has not been called due to a fault of the software or due to being truly absent. One main reason why standard variant callers
miss true variants is a low read coverage depth. Therefore we create a copy of the sequencing data with a smaller read coverage and perform variant calling on this
subsampled dataset. Next, variant calling on the original sequencing data (pre downsampling) is performed. Now, for each candidate variant there are the following options:

| called in subsampled data | called in original data | explanation |
|---------------------------|--------------------------|--------------|
| YES                     | YES                  | likely a true biological variant      |
| NO                     | NO                  | likely truly biologically absent      |
| NO                     | YES                  | likely biologically preset but missed by the variant caller      |
| YES                     | NO                  | something is very wrong      |

Therefore after running the pipeline on a large cohort of patients, the goal is to classify the "likely truly biologically absent" variants as TNs and the 
"likely biologically preset but missed by the variant caller" variants as FNs. The next step would be to train a classifier on such an annotated dataset.


## How to run
- create the conda environment based on ```save_the_environment.yaml```
- ensure that the ```snakedir/``` directory contains both the ```scripts/``` directory, the ```Snakefile.smk``` and the ```config.yaml``` file
- next, format your config file in the following way:
```
# data_________________________________________________________________________________________________________________
tumor_samples:    # paths to the multiple tumor samples
  sample_T1:    # the first tumor sample
    R1: T1_R1_001.fastq.gz  # the forward read
    R2: T1_R2_001.fastq.gz  # the reverse read
    tumor_sample_name: T1   # the sample name (needed for Mutect2)
  sample_T2:
    R1: T2_R1_001.fastq.gz
    R2: T2_R2_001.fastq.gz
    tumor_sample_name: T2

normal_sample:     # paths to the normal sample (only one normal sample per patient assumed)
  R1: N1_R1_001.fastq.gz
  R2: N1_R2_001.fastq.gz

normal_sample_name: N1

# parameters_________________________________________________________________________________________________________
snakedir: /snakedir      # the path to where your snakedir is
target_num_reads_tumor: 9719289    # number of reads to retain after subsampling in the tumor samples (more on that below)
target_num_reads_normal: 9719289   # number of reads to retain after subsampling in the normal sample  (more on that below)
subsampling_random_seed: 4      # the random seed for the random process of read subsampling

resources_mb:      # resources for each rule (in MB)
  subsample_tumor: 100000
  subsample_normal: 100000
  map_tumor_original: 100000
  map_tumor_subsampled: 100000
  map_normal_original: 100000
  map_normal_subsampled: 100000
  process_bam_tumor_original: 100000
  process_bam_tumor_subsampled: 100000
  process_bam_normal_original: 100000
  process_bam_normal_subsampled: 100000
  mutect2_call_original: 100000
  mutect2_call_subsampled: 100000
  mutect2_filter_original: 100000
  mutect2_filter_subsampled: 100000
  mutect2_keep_only_PASS_calls_original: 50000
  mutect2_keep_only_PASS_calls_subsampled: 50000
  collect: 50000
  multi_vcf_compilation: 25000

# genome
reference_genome: /snakedir/data/genome/hs37d5.fa  # the reference genome, note that all the index files for the referece must be present in the same directory
```

According to my knowledge there is no tool that subsampled reads up to a specified coverage. The tool used here ```seqtk``` subsampled up to a given number of reads.
How many reads you need to keep in order to attain your desired coverage is outlined in the file ```target_coverage_readme.txt```

Assuming your ```config.yaml``` file is correct, the conda environment is created and your ```snakedir/``` directory tree looks like this, you are ready to go.
```
├── Snakefile.smk
├── config.yaml
├── scripts
│   ├── MultiVcfVariantCompilation.py
│   ├── collect.py
│   ├── mapping.sh
│   ├── mutect2_call.sh
│   ├── mutect2_filter.sh
│   ├── mutect2_keep_only_PASS_calls.sh
│   ├── process_bams.sh
│   └── subsample_fastq.sh
└── target_coverage_readme.txt
```

#### Run
```
snakemake --snakefile Snakefile.smk --configfile config.yaml --cores <INT>
```

#### Run on BIH HPC
Here you can run this pipeline in the background as a slurm job. Make sure to run it on a login node and on a detachable screen (not literally).

```
snakemake --snakefile Snakefile.smk --configfile config.yaml --profile=cubi-v1 --jobs <INT>
```

After a succesfull run the directory tree of ```snakedir/``` might look like this:
```
├── Snakefile.smk
├── collected
│   └── collected_vcf_paths.txt
├── config.yaml
├── mapped_normal_original
│   └── N1.orig.bam
├── mapped_normal_subsampled
│   └── N1.subsampled.bam
├── mapped_processed_normal_original
│   ├── N1.orig.bam
│   ├── N1.orig.bam.bai
│   └── N1.orig.tmp.bam
├── mapped_processed_normal_subsampled
│   ├── N1.subsampled.bam
│   ├── N1.subsampled.bam.bai
│   └── N1.subsampled.tmp.bam
├── mapped_processed_tumor_original
│   ├── sample_T1.orig.bam
│   ├── sample_T1.orig.bam.bai
│   ├── sample_T1.orig.tmp.bam
│   ├── sample_T2.orig.bam
│   ├── sample_T2.orig.bam.bai
│   └── sample_T2.orig.tmp.bam
├── mapped_processed_tumor_subsampled
│   ├── sample_T1.subsampled.bam
│   ├── sample_T1.subsampled.bam.bai
│   ├── sample_T1.subsampled.tmp.bam
│   ├── sample_T2.subsampled.bam
│   ├── sample_T2.subsampled.bam.bai
│   └── sample_T2.subsampled.tmp.bam
├── mapped_tumor_original
│   ├── sample_T1.orig.bam
│   └── sample_T2.orig.bam
├── mapped_tumor_subsampled
│   ├── sample_T1.subsampled.bam
│   └── sample_T2.subsampled.bam
├── multi_vcf_compilation
│   └── multi_vcf_compilation.txt
├── mutect2_filtered_vcf_original
│   ├── sample_T1.orig.filtered.vcf.gz
│   ├── sample_T1.orig.filtered.vcf.gz.filteringStats.tsv
│   ├── sample_T1.orig.filtered.vcf.gz.tbi
│   ├── sample_T2.orig.filtered.vcf.gz
│   ├── sample_T2.orig.filtered.vcf.gz.filteringStats.tsv
│   └── sample_T2.orig.filtered.vcf.gz.tbi
├── mutect2_filtered_vcf_subsampled
│   ├── sample_T1.subsampled.filtered.vcf.gz
│   ├── sample_T1.subsampled.filtered.vcf.gz.filteringStats.tsv
│   ├── sample_T1.subsampled.filtered.vcf.gz.tbi
│   ├── sample_T2.subsampled.filtered.vcf.gz
│   ├── sample_T2.subsampled.filtered.vcf.gz.filteringStats.tsv
│   └── sample_T2.subsampled.filtered.vcf.gz.tbi
├── mutect2_only_PASS_vcf_original
│   ├── sample_T1.orig.PASS.vcf.gz
│   ├── sample_T1.orig.PASS.vcf.gz.tbi
│   ├── sample_T2.orig.PASS.vcf.gz
│   └── sample_T2.orig.PASS.vcf.gz.tbi
├── mutect2_only_PASS_vcf_subsampled
│   ├── sample_T1.subsampled.PASS.vcf.gz
│   ├── sample_T1.subsampled.PASS.vcf.gz.tbi
│   ├── sample_T2.subsampled.PASS.vcf.gz
│   └── sample_T2.subsampled.PASS.vcf.gz.tbi
├── mutect2_raw_vcf_original
│   ├── sample_T1.orig.vcf.gz
│   ├── sample_T1.orig.vcf.gz.stats
│   ├── sample_T1.orig.vcf.gz.tbi
│   ├── sample_T2.orig.vcf.gz
│   ├── sample_T2.orig.vcf.gz.stats
│   └── sample_T2.orig.vcf.gz.tbi
├── mutect2_raw_vcf_subsampled
│   ├── sample_T1.subsampled.vcf.gz
│   ├── sample_T1.subsampled.vcf.gz.stats
│   ├── sample_T1.subsampled.vcf.gz.tbi
│   ├── sample_T2.subsampled.vcf.gz
│   ├── sample_T2.subsampled.vcf.gz.stats
│   └── sample_T2.subsampled.vcf.gz.tbi
├── scripts
│   ├── MultiVcfVariantCompilation.py
│   ├── collect.py
│   ├── mapping.sh
│   ├── mutect2_call.sh
│   ├── mutect2_filter.sh
│   ├── mutect2_keep_only_PASS_calls.sh
│   ├── process_bams.sh
│   └── subsample_fastq.sh
├── subsampled_normal
│   ├── N1_R1.subsampled.fastq.gz
│   └── N1_R2.subsampled.fastq.gz
├── subsampled_tumor
│   ├── sample_T1_R1.subsampled.fastq.gz
│   ├── sample_T1_R2.subsampled.fastq.gz
│   ├── sample_T2_R1.subsampled.fastq.gz
│   └── sample_T2_R2.subsampled.fastq.gz
└── target_coverage_readme.txt
```

The final output file is ``` multi_vcf_compilation/ multi_vcf_compilation.txt``` and it has the following form:

|           | sample1.original | sample2.original | sample1.subsampled | sample2.subsampled |
|-----------|----------|----------|----------|----------|
| **locus1** | 1  | 1  | 1  | 1  |
| **locus2** | 1  | 1  | 0  | 1  |
| **locus3** | 0 | 1 | 0 | 1 |

In this example we could say that **locus2** is a FN (in sample1) and **locus3** is a TN (in sample1). 

