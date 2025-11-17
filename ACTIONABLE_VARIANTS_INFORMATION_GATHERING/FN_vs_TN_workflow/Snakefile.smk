import os
import re
 
# helper functions =======================================================================================================
def extract_basename(path):
    filename = os.path.basename(path)
    name = re.sub(r'\..*\.bam$', '', filename)
    return name

# ======================================================================================================================


# get data from the config file ======================================================================================
cf_snakedir = config['snakedir']
cf_reference_genome = config['reference_genome']
cf_tumor_sample_list = config['tumor_samples']
cf_normal_sample = config['normal_sample']
cf_normal_sample_name = config['normal_sample_name']
cf_target_num_reads_tumor = config['target_num_reads_tumor']
cf_target_num_reads_normal = config['target_num_reads_normal']
cf_random_subsampling_seed = config['subsampling_random_seed']
TUMORS = list(cf_tumor_sample_list.keys())  # list of tumor (fastq) samples (each sample has a forward and reverse read)
# below is a dict that maps each tumor samples to its sample name (e.g. sample_T1 --> T1)
cf_tumor_sample_names_dict = {sample:info['tumor_sample_name'] for sample, info in cf_tumor_sample_list.items()}
# =================================================================================================================


rule all:
    input:
        "multi_vcf_compilation/multi_vcf_compilation.txt"

rule subsample_tumor:
    resources:
        mem_mb=config['resources_mb']['subsample_tumor']
    input:
        R1 = lambda wc: cf_tumor_sample_list[wc.wildcard_tumor_sample]["R1"],
        R2 = lambda wc: cf_tumor_sample_list[wc.wildcard_tumor_sample]["R2"]
    output:
        R1 = "subsampled_tumor/{wildcard_tumor_sample}_R1.subsampled.fastq.gz",
        R2 = "subsampled_tumor/{wildcard_tumor_sample}_R2.subsampled.fastq.gz"
    params:
        rule_target_num_reads_tumor = cf_target_num_reads_tumor,
        rule_random_subsampling_seed = cf_random_subsampling_seed
    shell:
        """
        bash {cf_snakedir}/scripts/subsample_fastq.sh {params.rule_random_subsampling_seed} \
                                                    "{input.R1}" \
                                                    "{input.R2}" \
                                                    {params.rule_target_num_reads_tumor} \
                                                    "{output.R1}" \
                                                    {output.R2}
        """


rule subsample_normal:
    resources:
        mem_mb = config['resources_mb']['subsample_normal']
    input:
        R1 = cf_normal_sample["R1"],
        R2 = cf_normal_sample["R2"]
    output:
        R1 = "subsampled_normal/{cf_normal_sample_name}_R1.subsampled.fastq.gz",
        R2 = "subsampled_normal/{cf_normal_sample_name}_R2.subsampled.fastq.gz"
    params:
        rule_target_num_reads_normal = cf_target_num_reads_normal,
        rule_random_subsampling_seed = cf_random_subsampling_seed
    shell:
        """
        bash {cf_snakedir}/scripts/subsample_fastq.sh \
            {params.rule_random_subsampling_seed} \
            "{input.R1}" \
            "{input.R2}" \
            {params.rule_target_num_reads_normal} \
            "{output.R1}" \
            "{output.R2}"
        """



rule map_tumor_original:
    resources:
        mem_mb=config['resources_mb']['map_tumor_original']
    input:
        R1 = lambda wc: cf_tumor_sample_list[wc.wildcard_tumor_sample]["R1"],
        R2 = lambda wc: cf_tumor_sample_list[wc.wildcard_tumor_sample]["R2"]
    output:
        "mapped_tumor_original/{wildcard_tumor_sample}.orig.bam"
    params:
        threads=8,
        rule_reference_genome = cf_reference_genome
    shell:
         """
        bash {cf_snakedir}/scripts/mapping.sh \
            "{params.rule_reference_genome}" \
            "{input.R1}" \
            "{input.R2}" \
            "{output}" \
            {params.threads}
        """



rule map_normal_original:
    resources:
        mem_mb=config['resources_mb']['map_normal_original']
    input:
        R1 = cf_normal_sample["R1"],
        R2 = cf_normal_sample["R2"]
    output:
        "mapped_normal_original/{cf_normal_sample_name}.orig.bam"
    params:
        threads=8,
        rule_reference_genome = cf_reference_genome
    shell:
         """
        bash {cf_snakedir}/scripts/mapping.sh \
            "{params.rule_reference_genome}" \
            "{input.R1}" \
            "{input.R2}" \
            "{output}" \
            {params.threads}
        """



rule map_tumor_subsampled:
    resources:
        mem_mb=config['resources_mb']['map_tumor_subsampled']
    input:
        R1 = "subsampled_tumor/{wildcard_tumor_sample}_R1.subsampled.fastq.gz",
        R2 = "subsampled_tumor/{wildcard_tumor_sample}_R2.subsampled.fastq.gz"
    output:
        "mapped_tumor_subsampled/{wildcard_tumor_sample}.subsampled.bam"
    params:
        threads=8,
        rule_reference_genome = cf_reference_genome
    shell:
         """
        bash {cf_snakedir}/scripts/mapping.sh \
            "{params.rule_reference_genome}" \
            "{input.R1}" \
            "{input.R2}" \
            "{output}" \
            {params.threads}
        """


rule map_normal_subsampled:
    resources:
        mem_mb=config['resources_mb']['map_normal_subsampled']
    input:
        R1 = "subsampled_normal/{cf_normal_sample_name}_R1.subsampled.fastq.gz",
        R2 = "subsampled_normal/{cf_normal_sample_name}_R2.subsampled.fastq.gz"
    output:
        "mapped_normal_subsampled/{cf_normal_sample_name}.subsampled.bam"
    params:
        threads=8,
        rule_reference_genome = cf_reference_genome
    shell:
         """
        bash {cf_snakedir}/scripts/mapping.sh \
            "{params.rule_reference_genome}" \
            "{input.R1}" \
            "{input.R2}" \
            "{output}" \
            {params.threads}
        """


rule process_bam_tumor_original:
    resources:
        mem_mb=config['resources_mb']['process_bam_tumor_original']
    input:
        "mapped_tumor_original/{wildcard_tumor_sample}.orig.bam"
    output:
        output_bam = "mapped_processed_tumor_original/{wildcard_tumor_sample}.orig.bam",
        tmp_bam = "mapped_processed_tumor_original/{wildcard_tumor_sample}.orig.tmp.bam"
    params:
        threads=8,
        rule_tumor_sample_name=lambda wc: cf_tumor_sample_names_dict[wc.wildcard_tumor_sample]
    shell:
        """
        bash {cf_snakedir}/scripts/process_bams.sh "{input}" \
                                                    "{output.tmp_bam}" \
                                                    "{params.rule_tumor_sample_name}" \
                                                    "{params.threads}" \
                                                    "{output.output_bam}"
        """

rule process_bam_normal_original:
    resources:
        mem_mb=config['resources_mb']['process_bam_normal_original']
    input:
        "mapped_normal_original/{cf_normal_sample_name}.orig.bam"
    output:
        output_bam = "mapped_processed_normal_original/{cf_normal_sample_name}.orig.bam",
        tmp_bam = "mapped_processed_normal_original/{cf_normal_sample_name}.orig.tmp.bam"
    params:
        threads=8,
        rule_normal_sample_name=cf_normal_sample_name
    shell:
        """
        bash {cf_snakedir}/scripts/process_bams.sh "{input}" \
                                                    "{output.tmp_bam}" \
                                                    "{params.rule_normal_sample_name}" \
                                                    "{params.threads}" \
                                                    "{output.output_bam}"
        """

rule process_bam_tumor_subsampled:
    resources:
        mem_mb=config['resources_mb']['process_bam_tumor_subsampled']
    input:
        "mapped_tumor_subsampled/{wildcard_tumor_sample}.subsampled.bam"
    output:
        output_bam = "mapped_processed_tumor_subsampled/{wildcard_tumor_sample}.subsampled.bam",
        tmp_bam = "mapped_processed_tumor_subsampled/{wildcard_tumor_sample}.subsampled.tmp.bam"
    params:
        threads=8,
        rule_tumor_sample_name=lambda wc: cf_tumor_sample_names_dict[wc.wildcard_tumor_sample]
    shell:
        """
        bash {cf_snakedir}/scripts/process_bams.sh "{input}" \
                                                    "{output.tmp_bam}" \
                                                    "{params.rule_tumor_sample_name}" \
                                                    "{params.threads}" \
                                                    "{output.output_bam}"
        """

rule process_bam_normal_subsampled:
    resources:
        mem_mb=config['resources_mb']['process_bam_normal_subsampled']
    input:
        "mapped_normal_subsampled/{cf_normal_sample_name}.subsampled.bam"
    output:
        output_bam = "mapped_processed_normal_subsampled/{cf_normal_sample_name}.subsampled.bam",
        tmp_bam = "mapped_processed_normal_subsampled/{cf_normal_sample_name}.subsampled.tmp.bam"
    params:
        threads=8,
        rule_normal_sample_name=cf_normal_sample_name
    shell:
        """
        bash {cf_snakedir}/scripts/process_bams.sh "{input}" \
                                                    "{output.tmp_bam}" \
                                                    "{params.rule_normal_sample_name}" \
                                                    "{params.threads}" \
                                                    "{output.output_bam}"
        """


rule mutect2_call_original:
    resources:
        mem_mb=config['resources_mb']['mutect2_call_original']
    input:
        original_tumor_bam="mapped_processed_tumor_original/{wildcard_tumor_sample}.orig.bam",
        original_normal_bam=f"mapped_processed_normal_original/{cf_normal_sample_name}.orig.bam"
    output:
        "mutect2_raw_vcf_original/{wildcard_tumor_sample}.orig.vcf.gz"
    params:
        rule_reference_genome=cf_reference_genome,
        rule_normal_sample_name=cf_normal_sample_name,
        threads=8
    shell:
        """
        bash {cf_snakedir}/scripts/mutect2_call.sh "{input.original_tumor_bam}" \
                                                    "{input.original_normal_bam}" \
                                                    "{params.rule_normal_sample_name}" \
                                                    "{params.rule_reference_genome}" \
                                                    "{output}" \
                                                    "{params.threads}"
        """


rule mutect2_call_subsampled:
    resources:
        mem_mb=config['resources_mb']['mutect2_call_subsampled']
    input:
        subsampled_tumor_bam="mapped_processed_tumor_subsampled/{wildcard_tumor_sample}.subsampled.bam",
        subsampled_normal_bam=f"mapped_processed_normal_subsampled/{cf_normal_sample_name}.subsampled.bam"
    output:
        "mutect2_raw_vcf_subsampled/{wildcard_tumor_sample}.subsampled.vcf.gz"
    params:
        rule_reference_genome=cf_reference_genome,
        rule_normal_sample_name=cf_normal_sample_name,
        threads=8
    shell:
        """
        bash {cf_snakedir}/scripts/mutect2_call.sh "{input.subsampled_tumor_bam}" \
                                                    "{input.subsampled_normal_bam}" \
                                                    "{params.rule_normal_sample_name}" \
                                                    "{params.rule_reference_genome}" \
                                                    "{output}" \
                                                    "{params.threads}"
        """


rule mutect2_filter_original:
    resources:
        mem_mb=config['resources_mb']['mutect2_filter_original']
    input:
        "mutect2_raw_vcf_original/{wildcard_tumor_sample}.orig.vcf.gz"
    output:
        "mutect2_filtered_vcf_original/{wildcard_tumor_sample}.orig.filtered.vcf.gz"
    params:
        rule_reference_genome=cf_reference_genome
    shell:
        """
         bash {cf_snakedir}/scripts/mutect2_filter.sh "{params.rule_reference_genome}" \
                                                      "{input}" \
                                                      "{output}"
        """


rule mutect2_filter_subsampled:
    resources:
        mem_mb=config['resources_mb']['mutect2_filter_subsampled']
    input:
        "mutect2_raw_vcf_subsampled/{wildcard_tumor_sample}.subsampled.vcf.gz"
    output:
        "mutect2_filtered_vcf_subsampled/{wildcard_tumor_sample}.subsampled.filtered.vcf.gz"
    params:
        rule_reference_genome=cf_reference_genome
    shell:
        """
         bash {cf_snakedir}/scripts/mutect2_filter.sh "{params.rule_reference_genome}" \
                                                      "{input}" \
                                                      "{output}"
        """


rule mutect2_keep_only_PASS_calls_original:
    resources:
        mem_mb=config['resources_mb']['mutect2_keep_only_PASS_calls_original']
    input:
        "mutect2_filtered_vcf_original/{wildcard_tumor_sample}.orig.filtered.vcf.gz"
    output:
        "mutect2_only_PASS_vcf_original/{wildcard_tumor_sample}.orig.PASS.vcf.gz"
    shell:
        """
        bash {cf_snakedir}/scripts/mutect2_keep_only_PASS_calls.sh "{input}" \
                                                      "{output}"
        """



rule mutect2_keep_only_PASS_calls_subsampled:
    resources:
        mem_mb=config['resources_mb']['mutect2_keep_only_PASS_calls_subsampled']
    input:
        "mutect2_filtered_vcf_subsampled/{wildcard_tumor_sample}.subsampled.filtered.vcf.gz"
    output:
        "mutect2_only_PASS_vcf_subsampled/{wildcard_tumor_sample}.subsampled.PASS.vcf.gz"
    shell:
        """
        bash {cf_snakedir}/scripts/mutect2_keep_only_PASS_calls.sh "{input}" \
                                                      "{output}"
        """

rule collect:
    resources:
        mem_mb=config['resources_mb']['collect']
    input:
        original = expand("mutect2_only_PASS_vcf_original/{wildcard_tumor_sample}.orig.PASS.vcf.gz",wildcard_tumor_sample=TUMORS),
        subsampled = expand("mutect2_only_PASS_vcf_subsampled/{wildcard_tumor_sample}.subsampled.PASS.vcf.gz",wildcard_tumor_sample=TUMORS)
    output:
        "collected/collected_vcf_paths.txt"
    shell:
        """
        python3 {cf_snakedir}/scripts/collect.py {output} {input.original} {input.subsampled} 
        """

rule multi_vcf_compilation:
    resources:
        mem_mb=config['resources_mb']['multi_vcf_compilation']
    input:
        "collected/collected_vcf_paths.txt"
    output:
        "multi_vcf_compilation/multi_vcf_compilation.txt"
    shell:
        """
        python3 {cf_snakedir}/scripts/MultiVcfVariantCompilation.py --input {input} --output {output}
        """
    
