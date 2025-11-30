## Theoretical background

The goal of this PhD project is to increase the senistivity of variant calling of actionable variants. Therefore loci that are not being called by existing software are of interest.
If a given locus is not being called, it's due to these two main reasons:
  - it is truly biologically absent (no variant) --> true negative (TN)
  - it is a true variant that is missed by the software --> false negative (FN) --> **we want to capture these**

To achieve this a mathematical model is needed (most probably a Deep Learning model). Therefore training data, where the targets will be the TN and FN sites, is needed.
Since there exists no ground truth for this data, we need to guess whether a given locus that is not called is a TN or an FN. This is the task of the following workflow.

We have data from patients with multiple tumor biopsies.

<p align="center">
  <img src="https://github.com/JakubLiu/Supervised_detection_of_clinically_actionable_cancer_variants/blob/main/ACTIONABLE_VARIANTS_INFORMATION_GATHERING/FN_vs_TN_workflow_old/PIPELINE_V1/venn2.png.png" width="300">
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
| YES                     | NO                  | multiple reasons*      |

Therefore after running the pipeline on a large cohort of patients, the goal is to classify the "likely truly biologically absent" variants as TNs and the 
"likely biologically preset but missed by the variant caller" variants as FNs. The next step would be to train a classifier on such an annotated dataset.

*for example if (due to subsampling) some germline variants get missed in the normal, they can falsely appear as somatic variants in the tumor sample(s)

## How to run
- create the conda environment based on ```save_the_environment.yaml```
- ensure that the ```snakedir/``` directory contains both the ```scripts/``` directory, the ```Snakefile.smk``` and the ```config.yaml``` file
- next, format your config file in the following way:
```
# data_____________________________________________________________________________________________________
tumor_samples:
  sample_T1:
    R1: T1_R1_001.fastq.gz
    R2: T1_001.fastq.gz
    tumor_sample_name: T1
  sample_T2:
    R1: T2_R1_001.fastq.gz
    R2: T2_R2_001.fastq.gz
    tumor_sample_name: T2

normal_sample:
  R1: N1_R1_001.fastq.gz
  R2: N1_R2_001.fastq.gz
  normal_sample_name: N1

reference_genome: /data/genome/GRCh38.d1.vd1.fa
regions_file: /data/regions/Regions.bed


# parameters____________________________________________________________________________________________
snakedir: /snakedir
target_coverage_tumor: 30     
target_coverage_normal: 15   
subsampling_random_seed: 1
adapter: /data/adapters/TruSeq3-PE.fa

chromosomes:
  - chr1
  - chr2
  - chr3
  - chr4
  - chr5
  - chr6
  - chr7
  - chr8
  - chr9
  - chr10
  - chr11
  - chr12
  - chr13
  - chr14
  - chr15
  - chr16
  - chr17
  - chr18
  - chr19
  - chr20
  - chr21
  - chr22


# resources____________________________________________________________________________________________
resources_mb:
  adapter_trimming_tumor: 100000
  adapter_trimming_normal: 100000
  rule_map_tumor: 100000
  rule_map_normal: 100000
  rule_mark_duplicates_tumor: 100000
  rule_mark_duplicates_normal: 100000
  postprocess_bam_tumor: 100000
  postprocess_bam_normal: 100000
  samtools_depth_tumor: 100000
  samtools_depth_normal: 100000
  calculate_per_tumor_sample_depth: 50000
  calculate_subsampling_proportion_normal: 50000
  calculate_subsampling_proportion_tumor: 50000
  samtools_subsampled_tumor: 100000
  samtools_subsampled_normal: 100000
  process_subsampled_bam_tumor: 100000
  process_subsampled_bam_normal: 100000
  mutect2_scattercall_original: 100000
  mutect2_scattercall_subsampled: 100000
  filter_mutect2_calls_original: 100000
  filter_mutect2_calls_subsampled: 100000
  keep_PASS_original: 50000
  keep_PASS_subsampled: 50000
  collect: 25000
  multi_vcf_compilation: 100000
  gather_mutect2_stats_files_original: 50000
  gather_mutect2_stats_files_subsampled: 50000
  create_chrom_vcf_list_subsampled: 50000
  create_chrom_vcf_list_original: 50000
  bcftools_concat_original: 100000
  bcftools_concat_subsampled: 100000
  sort_index_original: 100000
  sort_index_subsampled: 100000
  mutect2_filter_original: 100000
  mutect2_filter_subsampled: 100000
  keep_PASS_original: 50000
  keep_PASS_subsampled: 50000
  collect: 25000
  multi_vcf_compilation: 100000


```

According to my knowledge there is no tool that subsampled reads up to a specified coverage. The tool used here ```seqtk``` subsampled up to a given number of reads.
How many reads you need to keep in order to attain your desired coverage is outlined in the file ```target_coverage_readme.txt```

Assuming your ```config.yaml``` file is correct, the conda environment is created and your ```snakedir/``` directory tree looks like this, you are ready to go.
```
├── Snakefile.smk
├── config.yaml
├── scripts
          ├── adapter_trimming.sh
          ├── bcftools_concat.sh
          ├── calc_subsampling_proportion_normal.py
          ├── calc_subsampling_proportion_tumor.py
          ├── collect.py
          ├── keep_PASS.sh
          ├── make_vcf_list.py
          ├── mapping.sh
          ├── mark_duplicates.sh
          ├── multi_vcf_compilation.py
          ├── mutect2_call.sh
          ├── mutect2_filter.sh
          ├── per_sample_depth.py
          ├── process_bams.sh
          ├── process_subsampled_bams.sh
          ├── samtools_depth.sh
          ├── samtools_subsample.sh
          ├── sort_index_vcf.sh
          ├── sum_stats_files.sh
          └── vcf_listing.sh
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

The final output file is ``` multi_vcf_compilation/ multi_vcf_compilation.txt``` and it has the following form:

|           | sample1.original | sample2.original | sample1.subsampled | sample2.subsampled |
|-----------|----------|----------|----------|----------|
| **locus1** | 1  | 1  | 1  | 1  |
| **locus2** | 1  | 1  | 0  | 1  |
| **locus3** | 0 | 1 | 0 | 1 |

In this example we could say that **locus2** is a FN (in sample1) and **locus3** is a TN (in sample1). 

