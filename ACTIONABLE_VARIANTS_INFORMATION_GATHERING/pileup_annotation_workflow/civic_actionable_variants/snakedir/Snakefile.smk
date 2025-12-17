cf_bamlist_file = config['bamfile_list']
cf_reference_genome = config['reference_genome']
cf_actionables_csv = config['actionables_csv']
cf_snpEff_dir = config['snpEff_dir']
cf_snakedir = config['snakedir']


rule all:
    input:
        'processed_pileup/pileup.processed.csv'

rule process_actionable_csv:
    resources:
        mem_mb=10000
    input:
        cf_actionables_csv
    output:
        vcf = 'actionables.vcf',
        bed = 'actionables.bed'
    shell:
        """
        bash {cf_snakedir}/scripts/csv_to_bed.sh "{input}" "{output.bed}"
        python3 {cf_snakedir}/scripts/csv_to_vcf.py "{input}" "{output.vcf}"
        """


rule annotate:
    resources:
        mem_mb=100000
    input:
        'actionables.vcf'
    output:
        'annotation/actionables.annotated.vcf'
    params:
        annotation_file = 'GRCh37.87',
        rule_snpEff_dir = cf_snpEff_dir
    shell:
        """
        bash {cf_snakedir}/scripts/annotation.sh "{input}" "{params.annotation_file}" "{output}" "{params.rule_snpEff_dir}"
        """


rule make_pileup:
    resources:
        mem_mb=10000
    input:
        cf_bamlist_file
    output:
        'pileup/pileup.txt'
    params:
        bed = 'actionables.bed',
        rule_reference_genome = cf_reference_genome
    shell:
        """
        bash {cf_snakedir}/scripts/pileup.sh "{input}" "{params.bed}" "{params.rule_reference_genome}" "{output}"
        """


rule process_pileup:
    resources:
        mem_mb = 100000
    input:
        actionables = cf_actionables_csv,
        pileup = 'pileup/pileup.txt'
    output:
        'processed_pileup/pileup.processed.csv'
    shell:
        """
        python3 {cf_snakedir}/scripts/process_pileup.py "{input.actionables}" "{input.pileup}" "{output}"
        """






