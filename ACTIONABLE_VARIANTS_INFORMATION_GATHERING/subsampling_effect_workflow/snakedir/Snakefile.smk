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
cf_normal_sample_name = config["normal_sample"]["normal_sample_name"]
cf_target_coverage_tumor = config['target_coverage_tumor']
cf_target_coverage_normal = config['target_coverage_normal']
cf_random_subsampling_seed = config['subsampling_random_seed']
cf_regions_bed = config['regions_file']   # in our case this is the exome capture kit bed file
cf_adapter = config['adapter']
TUMORS = list(cf_tumor_sample_list.keys())  # list of tumor (fastq) samples (each sample has a forward and reverse read)
CHROMS = config["chromosomes"]
# below is a dict that maps each tumor samples to its sample name (e.g. sample_T1 --> T1)
cf_tumor_sample_names_dict = {sample:info['tumor_sample_name'] for sample, info in cf_tumor_sample_list.items()}
# =================================================================================================================
 
rule all:
    input:
        "multi_vcf_compilation/multi_vcf_compilation.txt"
 



rule adapter_trimming_tumor:
    resources:
        mem_mb=config['resources_mb']['adapter_trimming_tumor']
    input:
        R1 = lambda wc: cf_tumor_sample_list[wc.wildcard_tumor_sample]["R1"],
        R2 = lambda wc: cf_tumor_sample_list[wc.wildcard_tumor_sample]["R2"]
    output:
        R1_trimmed_paired = 'trimmed_tumor/{wildcard_tumor_sample}_R1.trimmed.paired.fastq.gz',
        R2_trimmed_paired = 'trimmed_tumor/{wildcard_tumor_sample}_R2.trimmed.paired.fastq.gz',
        R1_trimmed_unpaired = 'trimmed_tumor/{wildcard_tumor_sample}_R1.trimmed.unpaired.fastq.gz',
        R2_trimmed_unpaired = 'trimmed_tumor/{wildcard_tumor_sample}_R2.trimmed.unpaired.fastq.gz'
    params:
        threads = 8,
        rule_adapter = cf_adapter
    shell:
        """
        bash {cf_snakedir}/scripts/adapter_trimming.sh \
            "{input.R1}" "{input.R2}" \
            "{output.R1_trimmed_paired}" "{output.R2_trimmed_paired}" \
            "{output.R1_trimmed_unpaired}" "{output.R2_trimmed_unpaired}" \
            {params.threads} \
            "{params.rule_adapter}"
        """


rule adapter_trimming_normal:
    resources:
        mem_mb=config['resources_mb']['adapter_trimming_normal']
    input:
        R1 = cf_normal_sample["R1"],
        R2 = cf_normal_sample["R2"]
    output:
        R1_trimmed_paired = 'trimmed_normal/{cf_normal_sample_name}_R1.trimmed.paired.fastq.gz',
        R2_trimmed_paired = 'trimmed_normal/{cf_normal_sample_name}_R2.trimmed.paired.fastq.gz',
        R1_trimmed_unpaired = 'trimmed_normal/{cf_normal_sample_name}_R1.trimmed.unpaired.fastq.gz',
        R2_trimmed_unpaired = 'trimmed_normal/{cf_normal_sample_name}_R2.trimmed.unpaired.fastq.gz'
    params:
        threads = 8,
        rule_adapter = cf_adapter
    shell:
        """
        bash {cf_snakedir}/scripts/adapter_trimming.sh \
            "{input.R1}" "{input.R2}" \
            "{output.R1_trimmed_paired}" "{output.R2_trimmed_paired}" \
            "{output.R1_trimmed_unpaired}" "{output.R2_trimmed_unpaired}" \
            {params.threads} \
            "{params.rule_adapter}"
        """


rule map_tumor:
    resources:
        mem_mb=config['resources_mb']['rule_map_tumor']
    input:
        R1 = 'trimmed_tumor/{wildcard_tumor_sample}_R1.trimmed.paired.fastq.gz',
        R2 = 'trimmed_tumor/{wildcard_tumor_sample}_R2.trimmed.paired.fastq.gz'
    output:
        'mapped_tumor/{wildcard_tumor_sample}.bam'
    params:
        threads=8,
        rule_reference_genome=cf_reference_genome
    shell:
        """
        bash {cf_snakedir}/scripts/mapping.sh \
            "{params.rule_reference_genome}" \
            "{input.R1}" \
            "{input.R2}" \
            "{output}" \
            {params.threads}
        """
 
 
rule map_normal:
    resources:
        mem_mb=config['resources_mb']['rule_map_normal']
    input:
        R1 = 'trimmed_normal/{cf_normal_sample_name}_R1.trimmed.paired.fastq.gz',
        R2 = 'trimmed_normal/{cf_normal_sample_name}_R2.trimmed.paired.fastq.gz'
    output:
        'mapped_normal/{cf_normal_sample_name}.bam'
    params:
        threads=8,
        rule_reference_genome=cf_reference_genome
    shell:
        """
        bash {cf_snakedir}/scripts/mapping.sh \
            "{params.rule_reference_genome}" \
            "{input.R1}" \
            "{input.R2}" \
            "{output}" \
            {params.threads}
        """
 
 
# this rule not only removes the duplicates but also sorts and indexes the bamfile
rule mark_duplicates_tumor:
    resources:
        mem_mb=config['resources_mb']['rule_mark_duplicates_tumor']
    input:
        'mapped_tumor/{wildcard_tumor_sample}.bam'
    output:
        'mapped_dedupa_tumor/{wildcard_tumor_sample}.dedupa.sorted.bam'
    params:
       
    shell:
        """
        bash {cf_snakedir}/scripts/mark_duplicates.sh \
                                        "{input}" \
                                        "{output}" \
                                        "{wildcards.wildcard_tumor_sample}"
        """
 
 
# this rule not only removes the duplicates but also sorts and indexes the bamfile
rule mark_duplicates_normal:
    resources:
        mem_mb=config['resources_mb']['rule_mark_duplicates_normal']
    input:
        'mapped_normal/{cf_normal_sample_name}.bam'
    output:
        'mapped_dedupa_normal/{cf_normal_sample_name}.dedupa.sorted.bam'
    shell:
        """
        bash {cf_snakedir}/scripts/mark_duplicates.sh \
                                        "{input}" \
                                        "{output}" \
                                        "{cf_normal_sample_name}"
        """
 
 
# add readgroups, sort, index
rule postprocess_bam_tumor:
    resources:
        mem_mb=config['resources_mb']['postprocess_bam_tumor']
    input:
        'mapped_dedupa_tumor/{wildcard_tumor_sample}.dedupa.sorted.bam'
    output:
        output_bam = "mapped_dedupa_processed_tumor/{wildcard_tumor_sample}.dedupa.rg.sorted.bam",
        tmp_bam = "mapped_dedupa_processed_tumor/{wildcard_tumor_sample}.dedupa.rg.sorted.tmp.bam"
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
 
 
# add readgroups, sort, index
rule postprocess_bam_normal:
    resources:
        mem_mb=config['resources_mb']['postprocess_bam_normal']
    input:
        'mapped_dedupa_normal/{cf_normal_sample_name}.dedupa.sorted.bam'
    output:
        output_bam = "mapped_dedupa_processed_normal/{cf_normal_sample_name}.dedupa.rg.sorted.bam",
        tmp_bam = "mapped_dedupa_processed_normal/{cf_normal_sample_name}.dedupa.rg.sorted.tmp.bam"
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
 
 
# just run samtools depth
rule samtools_depth_tumor:
    resources:
        mem_mb=config['resources_mb']['samtools_depth_tumor']
    input:
        "mapped_dedupa_processed_tumor/{wildcard_tumor_sample}.dedupa.rg.sorted.bam"
    output:
        "samtools_depth_tumor/{wildcard_tumor_sample}.depth.txt"
    params:
        rule_regions=cf_regions_bed
    shell:
        """
        bash {cf_snakedir}/scripts/samtools_depth.sh "{input}" "{output}" "{params.rule_regions}"
        """
 
# just run samtools depth
rule samtools_depth_normal:
    resources:
        mem_mb=config['resources_mb']['samtools_depth_normal']
    input:
        'mapped_dedupa_normal/{cf_normal_sample_name}.dedupa.sorted.bam'
    output:
        "samtools_depth_normal/{cf_normal_sample_name}.depth.txt"
    params:
        rule_regions=cf_regions_bed
    shell:
        """
        bash {cf_snakedir}/scripts/samtools_depth.sh "{input}" "{output}" "{params.rule_regions}"
        """
 
 
 
# calculate the median depth for each tumor sample
rule calculate_per_tumor_sample_depth:
    resources:
        mem_mb=config['resources_mb']['calculate_per_tumor_sample_depth']
    input:
        "samtools_depth_tumor/{wildcard_tumor_sample}.depth.txt"
    output:
        "per_tumor_sample_depth/{wildcard_tumor_sample}.txt"
    script:
        f"{cf_snakedir}/scripts/per_sample_depth.py"
 
 
# based on the depth and the target depth calculate the subsampling proportion
rule calculate_subsampling_proportion_normal:
    resources:
        mem_mb=config['resources_mb']['calculate_subsampling_proportion_normal']
    input:
        f"samtools_depth_normal/{cf_normal_sample_name}.depth.txt"
    output:
        "subsampling_proportion_normal/prop_normal.txt"
    params:
        rule_target_coverage_normal = cf_target_coverage_normal
    script:
        cf_snakedir + "/scripts/calc_subsampling_proportion_normal.py"
 
 
# based on the depth and the target depth calculate the subsampling proportion
rule calculate_subsampling_proportion_tumor:
    resources:
        mem_mb=config['resources_mb']['calculate_subsampling_proportion_tumor']
    input:
        expand("per_tumor_sample_depth/{tumor_sample}.txt", tumor_sample=TUMORS)
    output:
        "subsampling_proportion_tumor/prop_tumor.txt"
    params:
        rule_target_coverage_tumor = cf_target_coverage_tumor
    script:
        cf_snakedir + "/scripts/calc_subsampling_proportion_tumor.py"
 
 

rule samtools_subsampled_tumor:
    resources:
        mem_mb=config['resources_mb']['samtools_subsampled_tumor']
    input:
        bam="mapped_dedupa_processed_tumor/{wildcard_tumor_sample}.dedupa.rg.sorted.bam",
        subsampling_proportion_file="subsampling_proportion_tumor/prop_tumor.txt"
    output:
        "subsampled_tumor/{wildcard_tumor_sample}.subsampled.bam"
    params:
        rule_subsampling_seed = cf_random_subsampling_seed
    shell:
        """
        bash {cf_snakedir}/scripts/samtools_subsample.sh "{input.bam}" \
                                        "{input.subsampling_proportion_file}" \
                                        {params.rule_subsampling_seed} \
                                        "{output}"
        """
 
rule samtools_subsampled_normal:
    resources:
        mem_mb=config['resources_mb']['samtools_subsampled_normal']
    input:
        bam="mapped_dedupa_processed_normal/{cf_normal_sample_name}.dedupa.rg.sorted.bam",
        subsampling_proportion_file="subsampling_proportion_normal/prop_normal.txt"
    output:
        "subsampled_normal/{cf_normal_sample_name}.subsampled.bam"
    params:
        rule_subsampling_seed = cf_random_subsampling_seed
    shell:
        """
        bash {cf_snakedir}/scripts/samtools_subsample.sh "{input.bam}" \
                                        "{input.subsampling_proportion_file}" \
                                        {params.rule_subsampling_seed} \
                                        "{output}"
        """
 
 
rule process_subsampled_bam_tumor:
    resources:
        mem_mb=config['resources_mb']['process_subsampled_bam_tumor']
    input:
        "subsampled_tumor/{wildcard_tumor_sample}.subsampled.bam"
    output:
        "subsampled_tumor/{wildcard_tumor_sample}.subsampled.sorted.bam"
    shell:
        """
        bash {cf_snakedir}/scripts/process_subsampled_bams.sh "{input}" "{output}"
        """
 
 
rule process_subsampled_bam_normal:
    resources:
        mem_mb=config['resources_mb']['process_subsampled_bam_normal']
    input:
        "subsampled_normal/{cf_normal_sample_name}.subsampled.bam"
    output:
        "subsampled_normal/{cf_normal_sample_name}.subsampled.sorted.bam"
    shell:
        """
        bash {cf_snakedir}/scripts/process_subsampled_bams.sh "{input}" "{output}"
        """


rule mutect2_scattercall_subsampled:
    resources:
        mem_mb=config['resources_mb']['mutect2_scattercall_subsampled']
    input:
        tumor_bam = "subsampled_tumor/{wildcard_tumor_sample}.subsampled.sorted.bam",
        normal_bam =  f"subsampled_normal/{cf_normal_sample_name}.subsampled.sorted.bam",
    output:
        "mutect2_raw_vcf_subsampled/{wildcard_tumor_sample}_{chrom}.subsampled.vcf.gz"
    params:
        rule_reference_genome = cf_reference_genome,
        rule_normal_sample_name = cf_normal_sample_name,
        rule_chromosome = "{chrom}"
    shell:
        """
        bash {cf_snakedir}/scripts/mutect2_call.sh "{input.tumor_bam}" \
                                                    "{input.normal_bam}" \
                                                    "{params.rule_normal_sample_name}" \
                                                    "{params.rule_reference_genome}" \
                                                    "{params.rule_chromosome}" \
                                                    "{output}"
        """


rule mutect2_scattercall_original:
    resources:
        mem_mb=config['resources_mb']['mutect2_scattercall_original']
    input:
        tumor_bam = "mapped_dedupa_processed_tumor/{wildcard_tumor_sample}.dedupa.rg.sorted.bam",
        normal_bam =  f"mapped_dedupa_processed_normal/{cf_normal_sample_name}.dedupa.rg.sorted.bam",
    output:
        "mutect2_raw_vcf_original/{wildcard_tumor_sample}_{chrom}.original.vcf.gz"
    params:
        rule_reference_genome = cf_reference_genome,
        rule_normal_sample_name = cf_normal_sample_name,
        rule_chromosome = "{chrom}"
    shell:
        """
        bash {cf_snakedir}/scripts/mutect2_call.sh "{input.tumor_bam}" \
                                                    "{input.normal_bam}" \
                                                    "{params.rule_normal_sample_name}" \
                                                    "{params.rule_reference_genome}" \
                                                    "{params.rule_chromosome}" \
                                                    "{output}"
        """



rule gather_mutect2_stats_files_original:
    resources:
        mem_mb=config['resources_mb']['gather_mutect2_stats_files_original']
    input:
        lambda wc: expand("mutect2_raw_vcf_original/{sample}_{chrom}.original.vcf.gz", 
                          sample=wc.wildcard_tumor_sample, chrom=CHROMS)
    output:
        # I add the 'sorted' suffix, because it will be needed later
        "mutect2_gathered_original/{wildcard_tumor_sample}.concated.original.sorted.vcf.gz.stats"
    params:
        dirname="mutect2_raw_vcf_original",
        samplename=lambda wc: wc.wildcard_tumor_sample
    shell:
        """
        bash {cf_snakedir}/scripts/sum_stats_files.sh "{params.dirname}" "{params.samplename}" "{output}"
        """



rule gather_mutect2_stats_files_subsampled:
    resources:
        mem_mb=config['resources_mb']['gather_mutect2_stats_files_subsampled']
    input:
        lambda wc: expand("mutect2_raw_vcf_subsampled/{sample}_{chrom}.subsampled.vcf.gz", 
                          sample=wc.wildcard_tumor_sample, chrom=CHROMS)
    output:
        # I add the 'sorted' suffix because it will be needed later
        "mutect2_gathered_subsampled/{wildcard_tumor_sample}.concated.subsampled.sorted.vcf.gz.stats"
    params:
        dirname="mutect2_raw_vcf_subsampled",
        samplename=lambda wc: wc.wildcard_tumor_sample
    shell:
        """
        bash {cf_snakedir}/scripts/sum_stats_files.sh "{params.dirname}" "{params.samplename}" "{output}"
        """



rule create_chrom_vcf_list_original:
    resources:
        mem_mb=config['resources_mb']['create_chrom_vcf_list_original']
    input:
        "mutect2_gathered_original/{wildcard_tumor_sample}.concated.original.sorted.vcf.gz.stats"  # dummy inpute
    output:
        "mutect2_gathered_original/{wildcard_tumor_sample}.vcf_list.original.txt"
    params:
        dirname=os.path.abspath("mutect2_raw_vcf_original"),
        samplename=lambda wc: wc.wildcard_tumor_sample
    shell:
        """
        bash {cf_snakedir}/scripts/vcf_listing.sh "{params.dirname}" "{params.samplename}" "{output}"
        """


rule create_chrom_vcf_list_subsampled:
    resources:
        mem_mb=config['resources_mb']['create_chrom_vcf_list_subsampled']
    input:
        "mutect2_gathered_subsampled/{wildcard_tumor_sample}.concated.subsampled.sorted.vcf.gz.stats"  # dummy input
    output:
        "mutect2_gathered_subsampled/{wildcard_tumor_sample}.vcf_list.subsampled.txt"
    params:
        dirname=os.path.abspath("mutect2_raw_vcf_subsampled"),
        samplename=lambda wc: wc.wildcard_tumor_sample
    shell:
        """
        bash {cf_snakedir}/scripts/vcf_listing.sh "{params.dirname}" "{params.samplename}" "{output}"
        """



rule bcftools_concat_original:
    resources:
        mem_mb=config['resources_mb']['bcftools_concat_original']
    input:
        "mutect2_gathered_original/{wildcard_tumor_sample}.vcf_list.original.txt"
    output:
        "mutect2_gathered_original/{wildcard_tumor_sample}.concated.original.vcf.gz"
    shell:
        """
        bash {cf_snakedir}/scripts/bcftools_concat.sh "{input}" "{output}"
        """


rule bcftools_concat_subsampled:
    resources:
        mem_mb=config['resources_mb']['bcftools_concat_subsampled']
    input:
        "mutect2_gathered_subsampled/{wildcard_tumor_sample}.vcf_list.subsampled.txt"
    output:
        "mutect2_gathered_subsampled/{wildcard_tumor_sample}.concated.subsampled.vcf.gz"
    shell:
        """
        bash {cf_snakedir}/scripts/bcftools_concat.sh "{input}" "{output}"
        """


rule sort_index_original:
    resources:
        mem_mb=config['resources_mb']['sort_index_original']
    input:
        "mutect2_gathered_original/{wildcard_tumor_sample}.concated.original.vcf.gz"
    output:
        "mutect2_gathered_original/{wildcard_tumor_sample}.concated.original.sorted.vcf.gz"
    shell:
        """
        bash {cf_snakedir}/scripts/sort_index_vcf.sh "{input}" "{output}"
        """

rule sort_index_subsampled:
    resources:
        mem_mb=config['resources_mb']['sort_index_subsampled']
    input:
        "mutect2_gathered_subsampled/{wildcard_tumor_sample}.concated.subsampled.vcf.gz"
    output:
        "mutect2_gathered_subsampled/{wildcard_tumor_sample}.concated.subsampled.sorted.vcf.gz"
    shell:
        """
        bash {cf_snakedir}/scripts/sort_index_vcf.sh "{input}" "{output}"
        """


rule mutect2_filter_original:
    resources:
        mem_mb=config['resources_mb']['mutect2_filter_original']
    input:
        "mutect2_gathered_original/{wildcard_tumor_sample}.concated.original.sorted.vcf.gz"
    output:
        "mutect2_filtered_original/{wildcard_tumor_sample}.filtered.original.vcf.gz"
    params:
        rule_reference_genome=cf_reference_genome
    shell:
        """
        bash {cf_snakedir}/scripts/mutect2_filter.sh "{params.rule_reference_genome}" \
                                                    "{input}" "{output}"
        """

rule mutect2_filter_subsampled:
    resources:
        mem_mb=config['resources_mb']['mutect2_filter_subsampled']
    input:
        "mutect2_gathered_subsampled/{wildcard_tumor_sample}.concated.subsampled.sorted.vcf.gz"
    output:
        "mutect2_filtered_subsampled/{wildcard_tumor_sample}.filtered.subsampled.vcf.gz"
    params:
        rule_reference_genome=cf_reference_genome
    shell:
        """
        bash {cf_snakedir}/scripts/mutect2_filter.sh "{params.rule_reference_genome}" \
                                                    "{input}" "{output}"
        """


rule keep_PASS_original:
    resources:
        mem_mb=config['resources_mb']['keep_PASS_original']
    input:
        "mutect2_filtered_original/{wildcard_tumor_sample}.filtered.original.vcf.gz"
    output:
        "mutect2_PASS_original/{wildcard_tumor_sample}.PASS.original.vcf.gz"
    shell:
        """
        bash {cf_snakedir}/scripts/keep_PASS.sh "{input}" "{output}"
        """

rule keep_PASS_subsampled:
    resources:
        mem_mb=config['resources_mb']['keep_PASS_subsampled']
    input:
        "mutect2_filtered_subsampled/{wildcard_tumor_sample}.filtered.subsampled.vcf.gz"
    output:
        "mutect2_PASS_subsampled/{wildcard_tumor_sample}.PASS.subsampled.vcf.gz"
    shell:
        """
        bash {cf_snakedir}/scripts/keep_PASS.sh "{input}" "{output}"
        """


rule collect:
    resources:
        mem_mb=config['resources_mb']['collect']
    input:
        original = expand("mutect2_PASS_original/{wildcard_tumor_sample}.PASS.original.vcf.gz", wildcard_tumor_sample=TUMORS),
        subsampled = expand("mutect2_PASS_subsampled/{wildcard_tumor_sample}.PASS.subsampled.vcf.gz", wildcard_tumor_sample = TUMORS)
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
        python3 {cf_snakedir}/scripts/multi_vcf_compilation.py --input {input} --output {output}
        """


