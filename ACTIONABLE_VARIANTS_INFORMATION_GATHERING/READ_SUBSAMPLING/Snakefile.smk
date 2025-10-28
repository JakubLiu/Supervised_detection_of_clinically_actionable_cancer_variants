import os

def extract_basename(path):
    return os.path.splitext(os.path.basename(path))[0]


bamfiles = config['bam_files']
subsample_proportions = config['proportion_reads_to_keep']
ref = config['reference_genome']

# this rule (the expand command) defines all the (file x proportion) combinations that need to be created by the other rules
rule all:
    input:
        expand(
            "filtered_mutect2_calls/{bam_basename}_proportion_reads_kept_{prop}.filtered.vcf.gz",
            bam_basename=[extract_basename(f) for f in bamfiles],
            prop=subsample_proportions
        )


rule subsample_sort_index:
    output:
        bam = "subsampled_reads/{bam_basename}_proportion_reads_kept_{prop}.bam",
        bai = "subsampled_reads/{bam_basename}_proportion_reads_kept_{prop}.bam.bai"

    input:
        lambda wildcards: next(
            f for f in bamfiles if extract_basename(f) == wildcards.bam_basename
        )

    params:
        prop = "{prop}"

    shell:
        """
        samtools view -s {params.prop} -b {input} | samtools sort -o {output.bam}
        samtools index {output.bam}
        """


rule mutect2_call:
    input:
        bam = "subsampled_reads/{bam_basename}_proportion_reads_kept_{prop}.bam",
        bai = "subsampled_reads/{bam_basename}_proportion_reads_kept_{prop}.bam.bai"

    output:
        vcf = "raw_mutect2_calls/{bam_basename}_proportion_reads_kept_{prop}.unfiltered.vcf.gz"

    params:
        threads = 8,
        ref = ref

    shell:
        """
        gatk Mutect2 \
            -R {params.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            --native-pair-hmm-threads {params.threads}
        """


rule mutect2_filter:
    input:
        vcf = "raw_mutect2_calls/{bam_basename}_proportion_reads_kept_{prop}.unfiltered.vcf.gz",
        ref = ref

    output:
        vcf = "filtered_mutect2_calls/{bam_basename}_proportion_reads_kept_{prop}.filtered.vcf.gz"

    shell:
        """
        gatk FilterMutectCalls -R {input.ref} -V {input.vcf} -O {output.vcf}
        """

