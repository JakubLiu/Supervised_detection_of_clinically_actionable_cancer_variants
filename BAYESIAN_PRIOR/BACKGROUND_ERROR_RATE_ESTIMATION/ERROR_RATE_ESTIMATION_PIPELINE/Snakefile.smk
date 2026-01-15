import os

def out_bam(bam):
    return os.path.join(
        "masked",
        os.path.basename(bam).replace(".bam", ".masked.bam")
    )


cf_bamlist_file = config['bamlist_file']
cf_reference_genome = config['reference_genome']
cf_snakedir = config['snakedir']
cf_bamfile_extension = config['bamfile_extension']
cf_gene_bedfile = config['gene_bedfile']
cf_gene_name = config['gene_name']
cf_min_baseQ = config['minQ']

with open(cf_bamlist_file, 'r') as f:
    BAMS = [line.strip() for line in f if line.strip()]

SAMPLES = [
    os.path.basename(bam).replace(".bam", "")
    for bam in BAMS
]


BAMDICT = {
    os.path.basename(bam).replace(".bam", ""): bam
    for bam in BAMS
}


rule all:
    input:
        f'error_rate/{cf_gene_name}.error_rates.txt',
        expand(
            "masked/{sample}.masked.bam",
            sample=SAMPLES
        )



rule mask_lowQ_bases:
    resources:
        mem_mb=config['resources_mb']['rule_mask_lowQ_bases']
    input:
        bam=lambda wildcards: BAMDICT[wildcards.sample]
    output:
        'masked/{sample}.masked.bam'
    params:
        minQ = cf_min_baseQ
    shell:
        """
        bash {cf_snakedir}/scripts/mask_bases.exe "{input.bam}" "{output}" "{params.minQ}"
        """


rule create_masked_bamlist:
    resources:
        mem_mb=config['resources_mb']['rule_create_masked_bamlist']
    input:
        expand("masked/{sample}.masked.bam", sample=SAMPLES)
    output:
        'masked/masked_bamlist.txt'
    shell:
        """
        pritnf "%s\n" {input} > {output}
        """

rule pileup:
    resources:
        mem_mb = config['resources_mb']['rule_pileup']
    input:
        'masked/masked_bamlist.txt'
    output:
        f'pileup/{cf_gene_name}.pileup.txt'
    params:
        gene_bedfile = cf_gene_bedfile,
        reference_genome = cf_reference_genome
    shell:
         """
        bash {cf_snakedir}/scripts/pileup.sh "{input}" "{params.gene_bedfile}" "{params.reference_genome}" "{output}"
        rm masked/masked_bamlist.txt
        """



rule process_pileup:
    resources:
        mem_mb = config['resources_mb']['rule_process_pileup']
    input:
        f'pileup/{cf_gene_name}.pileup.txt'
    output:
        f'procesed_pileup/{cf_gene_name}.pileup.processed.txt'
    params:
        bamlist_file = cf_bamlist_file,
        bamfile_extension = cf_bamfile_extension
    shell:
        """
        python3 {cf_snakedir}/scripts/process_pileup.py "{input}" "{output}" "{params.bamlist_file}" "{params.bamfile_extension}"
        """

rule estimate_error_rate:
    resources:
        mem_mb = config['resources_mb']['rule_estimate_error_rate']
    input:
        f'procesed_pileup/{cf_gene_name}.pileup.processed.txt'
    output:
        f'error_rate/{cf_gene_name}.error_rates.txt'
    shell:
        """
        python3 {cf_snakedir}/scripts/error_rate.py "{input}" "{output}"
        """