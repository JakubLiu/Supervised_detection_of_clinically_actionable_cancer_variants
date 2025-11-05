import os
 
def extract_basename(path):
    return os.path.splitext(os.path.basename(path))[0]

# based on the --config remove_intermediate_dirs remove or keep the intermediate directories
def maybe_temp(path):
    if remove_intermediate == True:
        return temp(path)
    else:
        return path


# the directory from which snakemake has been run (otherwise there are errors when running it on an HPC)
rundir = config["source_dir"]


# by default remove all intermadiate directories
# can be changed via snakemake --config remove_intermediate_dirs=False [... other params] 
remove_intermediate = config.get("remove_intermediate_dirs", True)

reference_genome = config['reference_genome']   # make sure all the supporting files for the reference are there
bamfile_list = config['bamfile_list']     # list of bamfiles, one file per line
cutoffs = config['cutoffs']              # for example classify a locus as a FN when >= 10% of the reads support the alternative
coverages = config['coverages']         # the proportion of reads to keep durign bam file subsampling
sample_proportions = config['sample_proportions']   # for example fetch variants that where missed in 1/5 samples
actionable_variants_bed = config['actionable_variants_bed']  # bed file with actionable variants (e.g. from Civic)
 
 


rule all:     # modify this rule as you go on
    input:
      expand('postprocessed_mpileup/{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}_in_{sample_prop}_samples_ACTIONALBES_postprocessed_pileup.txt',
      raw_bam = [extract_basename(f) for f in bamfile_list],
      prop_reads_to_keep = coverages,
      sample_prop = sample_proportions)
 


rule subsample_bams:
    resources:
        mem_mb=32000
    
    input:
        lambda wildcards: next(f for f in bamfile_list if extract_basename(f)==wildcards.raw_bam)
 
    output:
        bam = maybe_temp('subsampled_bams/{raw_bam}_cov_{prop_reads_to_keep}.subsampled.bam'),
        bai = maybe_temp('subsampled_bams/{raw_bam}_cov_{prop_reads_to_keep}.subsampled.bam.bai')
 
    params:
        prop_reads_to_keep='{prop_reads_to_keep}'
 
    shell:
        '''
        samtools view -s {params.prop_reads_to_keep} -b {input} | samtools sort -o {output.bam}
        samtools index {output.bam}
        '''



rule mutect2_call:
    resources:
        mem_mb=32000
    input:
        bam = 'subsampled_bams/{raw_bam}_cov_{prop_reads_to_keep}.subsampled.bam',
        bai = 'subsampled_bams/{raw_bam}_cov_{prop_reads_to_keep}.subsampled.bam.bai',
        ref=reference_genome
    output:
        vcf = maybe_temp("raw_mutect2_calls/{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}.unfiltered.vcf.gz")
    params:
        threads=8
    shell:
        f"bash {rundir}/scripts/mutect2_call.sh {{input.bam}} {{input.ref}} {{output.vcf}} {{params.threads}}"




rule mutect2_filter:
    resources:
        mem_mb=32000
    input:
        vcf = "raw_mutect2_calls/{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}.unfiltered.vcf.gz",
        ref=reference_genome
    output:
        maybe_temp("filtered_mutect2_calls/{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}.vcf.gz")
    shell:
        f'bash {rundir}/scripts/mutect2_filter.sh {{input.ref}} {{input.vcf}} {{output}}'



rule keep_PASS_variant_calls:
    resources:
        mem_mb=3200
    input:
        "filtered_mutect2_calls/{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}.vcf.gz"
    output:
        maybe_temp(os.path.join(rundir, "only_PASS_variant_calls", "{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}.onlyPASS.vcf.gz"))
    shell:
        f'bash {rundir}/scripts/keep_only_PASS_variant_calls.sh {{input}} {{output}}'




rule compile_multiple_vcfs:
    resources:
        mem_mb=32000
    input:
        expand(
            os.path.join(rundir, "only_PASS_variant_calls", "{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}.onlyPASS.vcf.gz"),
            raw_bam=[extract_basename(f) for f in bamfile_list],
            prop_reads_to_keep=coverages
        )
    output:
        maybe_temp(os.path.join(rundir, "multi_vcf_compilation", "multi_vcf_compilation.csv"))
    shell:
        f"""
        ls -l {rundir}/only_PASS_variant_calls/*.onlyPASS.vcf.gz | awk '{{{{print $NF}}}}' > vcf_filelist.txt
        python3 {rundir}/scripts/MultiVcfVariantCompilation.py --input vcf_filelist.txt --output {{output}}
        """


rule n_out_of_m_samples:
    resources:
        mem_mb=32000
    input:
        os.path.join(rundir, "multi_vcf_compilation", "multi_vcf_compilation.csv")
    params:
        sample_prop="{sample_prop}"
    output:
        maybe_temp("n_out_of_m_samples/{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}_in_{sample_prop}_samples.bed")
    shell:
        f'python3 {rundir}/scripts/VariantExtractor.py --input {{input}} --output {{output}} --proportion {{params.sample_prop}}'

rule merge_with_actionables_bed:
    resources:
        mem_mb=3200
    input:
        n_out_of_m_bed="n_out_of_m_samples/{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}_in_{sample_prop}_samples.bed",
        actionables_bed=actionable_variants_bed
    output:
        maybe_temp("merged_bed/{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}_in_{sample_prop}_samples_ACTIONALBES.bed")
    shell:
        """
        cat {input.n_out_of_m_bed} {input.actionables_bed} > {output}
        """


# this returns one mpileup per bamfile (not one mpileup per multiple bamfiles)
rule mpileup:
    resources:
        mem_mb=3200
    input:
        merged_bed="merged_bed/{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}_in_{sample_prop}_samples_ACTIONALBES.bed",
        bam='subsampled_bams/{raw_bam}_cov_{prop_reads_to_keep}.subsampled.bam',
        bai='subsampled_bams/{raw_bam}_cov_{prop_reads_to_keep}.subsampled.bam.bai',
        ref=reference_genome
    output:
        maybe_temp('samtools_mpileup/{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}_in_{sample_prop}_samples_ACTIONALBES_raw_pileup.txt')
    shell:
        """
        samtools mpileup -f {input.ref} {input.bam} -l {input.merged_bed} -s > {output}
        """

rule postprocess_mpileup:
    resources:
        mem_mb=3200
    input:
        'samtools_mpileup/{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}_in_{sample_prop}_samples_ACTIONALBES_raw_pileup.txt',
    output:
        maybe_temp('postprocessed_mpileup/{raw_bam}_proportion_reads_kept_{prop_reads_to_keep}_in_{sample_prop}_samples_ACTIONALBES_postprocessed_pileup.txt')
    script:
        f"{rundir}/scripts/postprocess_mpileup_FN_call.py"
 
