import os

cf_reference_genome = config['reference_genome']
cf_gene_bed = config['gene_bed']
cf_bamlist_file = config['bamlist_file']
cf_gzvcf_list_file = config['gzvcf_list_file']
cf_snakedir = config['snakedir']

with open(cf_bamlist_file) as f:
    BAM_FILES = [line.strip() for line in f if line.strip()]

with open(cf_gzvcf_list_file) as f:
    GZVCF_FILES =  [line.strip() for line in f if line.strip()]

VCF_TO_BAM_DICT = {}


if len(BAM_FILES) != len(GZVCF_FILES):
    raise ValueError(
        f"VCF list ({len(GZVCF_FILES)}) and BAM list ({len(BAM_FILES)}) have different lengths"
    )

for vcf, bam in zip(GZVCF_FILES, BAM_FILES):
    VCF_TO_BAM_DICT[vcf] = bam

# this dict frees me from having to use two different IDs for the bams and vcfs
PAIR_IDS = [os.path.basename(vcf) for vcf in VCF_TO_BAM_DICT] 

wildcard_constraints:
    vcf=".+"


# ------------------------------------------------------------------------------------------------------------------------------------------------------------


rule all:
    input:
        expand(
            "pileup_missed_loci_processed/{pair_id}.pileup.processed.txt",
            pair_id=PAIR_IDS
        )


rule identify_missing_loci:
    resources:
        mem_mb = config['resources']['rule_identify_missing_loci']
    input:
        rule_gene_bed = cf_gene_bed,
        rule_gzvcf = lambda wc: next(
            v for v in VCF_TO_BAM_DICT if os.path.basename(v) == wc.pair_id
        )
    output:
        'missed_loci/{pair_id}.missed_loci.txt'
    shell:
        """
        bash {cf_snakedir}/scripts/IdentifyMissedLoci.sh  "{input.rule_gene_bed}" "{input.rule_gzvcf}" "{output}"
        """


# this can take only one bamfile at a time
rule make_bam_readcount:
    resources:
        mem_mb = config['resources']['rule_make_bam_readcount']
    input:
        bam = lambda wc: VCF_TO_BAM_DICT[
            next(v for v in VCF_TO_BAM_DICT if os.path.basename(v) == wc.pair_id)
        ],
        bed = 'missed_loci/{pair_id}.missed_loci.txt'
    output:
        'pileup_missed_loci_raw/{pair_id}.pileup.raw.txt'
    params:
        rule_reference_genome = cf_reference_genome
    shell:
        """
        bash {cf_snakedir}/scripts/make_bam_readcounts.sh "{params.rule_reference_genome}" "{input.bed}" "{input.bam}" "{output}"
        """


rule process_bam_readcount:
    resources:
        mem_mb = config['resources']['rule_process_bam_readcount']
    input:
        'pileup_missed_loci_raw/{pair_id}.pileup.raw.txt'
    output:
        'pileup_missed_loci_processed/{pair_id}.pileup.processed.txt'
    shell:
        """
        python3 {cf_snakedir}/scripts/process_bam_readcount.py "{input}" "{output}"
        """