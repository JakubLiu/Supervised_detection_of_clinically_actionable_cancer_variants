# Warning
**This pipeline produces valid results only for SNV's, for indels or other variants it is not stable!**
# How to run
## required input files
**multisample calling file**
This file shows whether a given locus has been called (as a variant) in the given tumor sample or not.
(more on this file format here: https://github.com/JakubLiu/Supervised_detection_of_clinically_actionable_cancer_variants/tree/main/ACTIONABLE_VARIANTS_INFORMATION_GATHERING/subsampling_effect_workflow/snakedir)
Note that this pipeline requires the multisample calling file to have exactly one column per tumor sample and a locus column.
Example:
| locus | sample T1 | sample T2 | sample T3 |
|----------|----------|----------|----------|
| chr;pos;ref;alt  | 1 | 1 | 1 |
| chr;pos;ref;alt | 0 | 0 | 1 |
| chr;pos;ref;alt | 1 | 1 | 0 |
| chr;pos;ref;alt | 1 | 0 | 1 |

**bamfile list**
Just a text file with paths to the bamfiles of the tumor samples (one path per line). Note that the order of the bamfiles must
be identical to the order of the tumor samples in the multisample calling file

Some input files (and some other parameters) must be specified must be specified in the ```config.yaml```.

```
# parameters____________________________________________________________________________________________
snakedir: /snakedir
snpEff_dir: /snpeff/snpEff
reference_genome: /genome/GRCh38.d1.vd1.fa
snpEff_annotation_file: GRCh38.92
bamlist_file: /snakedir/init_input_files/bamlist.txt
 
 
# resources____________________________________________________________________________________________
resources_mb:
  make_vcf_to_annotate: 50000
  annotate: 100000
  pileup: 100000
  process_pileup: 50000
  mutate_and_merge: 50000

```

## Other requirements
**snpEff** installation with the correct annotation file.
For more information check: https://github.com/JakubLiu/Supervised_detection_of_clinically_actionable_cancer_variants/tree/main/ACTIONABLE_VARIANTS_INFORMATION_GATHERING/pileup_annotation_workflow/civic_actionable_variants
or: https://pcingola.github.io/SnpEff/

**Conda environment**
The environemnt can be created by using the ```environmentalists.yaml``` file.

**Before running** the pipeline, the directory tree should look like this:

```
├── Snakefile.smk
├── config.yaml
├── init_input_files
│   ├── bamlist.txt
│   └── multi_vcf_compilation.txt
├── scripts
   ├── annotate.sh
   ├── extract_loci.sh
   ├── make_vcf.py
   ├── mutate_and_merge.py
   ├── pileup.sh
   └── process_pileup.py

```

**After running** the pipeline, the directory tree should look like this:

```
├── Snakefile.smk
├── annotation
│   ├── annotated.vcf
│   └── vcf_to_annotate.vcf
├── config.yaml
├── environmentalists.yaml
├── init_input_files
│   ├── bamlist.txt
│   └── multi_vcf_compilation.txt
├── output
│   └── final_output.csv
├── pileup
│   ├── pileup.processed.csv
│   └── raw_pileup.txt
├── scripts
│   ├── annotate.sh
│   ├── extract_loci.sh
│   ├── make_vcf.py
│   ├── mutate_and_merge.py
│   ├── pileup.sh
│   └── process_pileup.py
├── snpEff_genes.txt
├── snpEff_summary.html
└── variants_loci
    └── loci.txt

```

## Command to run
```
snakemake --snakefile Snakefile.smk --configfile config.yaml --cores <num cores>
```

# Output files
## Annotation file
Just a VCF file of the variants from the multisample calling file, annotated using snpEff.

## Enchanced multisample calling file
| chrom | pos | ref | variant_alt | sample T1 | sample T2 | coverage sample T1 | coverage sample T2 | total alt reads sample T1 | total alt reads sample T2 | variant supporting alt reads sample T1 | variant supporting alt reads sample T2 | pileup sample T1 | pileup sample T2 |
|------|------|------|------|------|------|------|------|------|-------|-------|-------|-------|-------|
| 1 | 12 | A | G | 1 | 1 | 32 | 26 | 4 | 3 | 3 | 3 | AAAA.... | AAAGA... |
| 2 | 112 | A | T | 0 | 1 | 56 | 45 | 4 | 8 | 4 | 7 | ATAA.... | AAAAA... |
