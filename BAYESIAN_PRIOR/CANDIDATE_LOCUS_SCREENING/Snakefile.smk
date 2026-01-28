import os

cf_reference_genome = config['reference_genome']
cf_gene_bed = config['gene_bed']
cf_bamlist_file_tumor = config['bamlist_file_tumor']
cf_bamlist_file_normal = config['bamlist_file_normal']
cf_gzvcf_list_file = config['gzvcf_list_file']
cf_snakedir = config['snakedir']

with open(cf_bamlist_file_tumor) as f:
    BAM_FILES_TUMOR = [line.strip() for line in f if line.strip()]


with open(cf_bamlist_file_normal) as f:
    BAM_FILES_NORMAL = [line.strip() for line in f if line.strip()]

with open(cf_gzvcf_list_file) as f:
    GZVCF_FILES =  [line.strip() for line in f if line.strip()]

VCF_TO_BAM_DICT_TUMOR = {}
if len(BAM_FILES_TUMOR) != len(GZVCF_FILES):
    raise ValueError(
        f"VCF list ({len(GZVCF_FILES)}) and BAM list ({len(BAM_FILES_TUMOR)}) have different lengths"
    )
for vcf, bam in zip(GZVCF_FILES, BAM_FILES_TUMOR):
    VCF_TO_BAM_DICT_TUMOR[vcf] = bam



VCF_TO_BAM_DICT_NORMAL = {}
if len(BAM_FILES_NORMAL) != len(GZVCF_FILES):
    raise ValueError(
        f"VCF list ({len(GZVCF_FILES)}) and BAM list ({len(BAM_FILES_NORMAL)}) have different lengths"
    )
for vcf, bam in zip(GZVCF_FILES, BAM_FILES_NORMAL):
    VCF_TO_BAM_DICT_NORMAL[vcf] = bam

# this dict frees me from having to use two different IDs for the bams and vcfs
PAIR_IDS = [os.path.basename(vcf) for vcf in VCF_TO_BAM_DICT_TUMOR] 

wildcard_constraints:
    vcf=".+"


# ------------------------------------------------------------------------------------------------------------------------------------------------------------


rule all:
    input:
        expand(
            'merged/{pair_id}.merged.csv',
            pair_id=PAIR_IDS
        )

rule identify_missing_loci:
    resources:
        mem_mb = config['resources']['rule_identify_missing_loci']
    input:
        rule_gene_bed = cf_gene_bed,
        rule_gzvcf = lambda wc: next(
            v for v in VCF_TO_BAM_DICT_TUMOR if os.path.basename(v) == wc.pair_id
        )
    output:
        'missed_loci/{pair_id}.missed_loci.txt'
    shell:
        """
        bash {cf_snakedir}/scripts/IdentifyMissedLoci.sh  "{input.rule_gene_bed}" "{input.rule_gzvcf}" "{output}"
        """


# make a read pileup in the tumor sample at the location of the missed loci
rule make_bam_readcount_tumor:
    resources:
        mem_mb = config['resources']['rule_make_bam_readcount_tumor']
    input:
        bam = lambda wc: VCF_TO_BAM_DICT_TUMOR[
            next(v for v in VCF_TO_BAM_DICT_TUMOR if os.path.basename(v) == wc.pair_id)
        ],
        bed = 'missed_loci_tumor/{pair_id}.missed_loci.tumor.txt'
    output:
        'pileup_missed_loci_raw_tumor/{pair_id}.pileup.raw.tumor.txt'
    params:
        rule_reference_genome = cf_reference_genome
    shell:
        """
        bash {cf_snakedir}/scripts/make_bam_readcounts.sh "{params.rule_reference_genome}" "{input.bed}" "{input.bam}" "{output}"
        """


# make a read pileup in the normal sample at the location of the missed loci
rule make_bam_readcount_normal:
    resources:
        mem_mb = config['resources']['rule_make_bam_readcount_normal']
    input:
        bam = lambda wc: VCF_TO_BAM_DICT_NORMAL[
            next(v for v in VCF_TO_BAM_DICT_NORMAL if os.path.basename(v) == wc.pair_id)
        ],
        bed = 'missed_loci/{pair_id}.missed_loci.txt'
    output:
        'pileup_missed_loci_raw_normal/{pair_id}.pileup.raw.normal.txt'
    params:
        rule_reference_genome = cf_reference_genome
    shell:
        """
        bash {cf_snakedir}/scripts/make_bam_readcounts.sh "{params.rule_reference_genome}" "{input.bed}" "{input.bam}" "{output}"
        """

# process the tumor pileup
rule process_bam_readcount_tumor:
    resources:
        mem_mb = config['resources']['rule_process_bam_readcount_tumor']
    input:
        'pileup_missed_loci_raw/{pair_id}.pileup.raw.txt'
    output:
        'pileup_missed_loci_processed/{pair_id}.pileup.processed.txt'
    shell:
        """
        python3 {cf_snakedir}/scripts/process_bam_readcount.py "{input}" "{output}"
        """



# process the tumor pileup
rule process_bam_readcount_normal:
    resources:
        mem_mb = config['resources']['rule_process_bam_readcount_normal']
    input:
        'pileup_missed_loci_raw_normal/{pair_id}.pileup.raw.normal.txt'
    output:
        'pileup_missed_loci_processed_normal/{pair_id}.pileup.processed.normal.txt'
    shell:
        """
        python3 {cf_snakedir}/scripts/process_bam_readcount.py "{input}" "{output}"
        """


rule merge:
    resources:
        mem_mb = config['resources']['rule_merge']
    input:
        tumor = 'pileup_missed_loci_processed/{pair_id}.pileup.processed.txt',
        normal = 'pileup_missed_loci_processed_normal/{pair_id}.pileup.processed.normal.txt'
    output:
        'merged/{pair_id}.merged.csv'
    shell:
        """
        python3 {cf_snakedir}/scripts/merge.py "{input.tumor}" "{input.normal}" "{output}"
        """