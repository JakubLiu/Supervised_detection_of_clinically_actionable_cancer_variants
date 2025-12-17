cf_bamlist_file = config['bamfile_list']
cf_reference_genome = config['reference_genome']
cf_actionables_csv = config['actionables_csv']
cf_snpEff_dir = config['snpEff_dir']
cf_snakedir = config['snakedir']


rule all:
    input:
        'processed_pileup/pileup.processed.csv',
        'annotation/actionables.annotated.vcf'

rule process_actionable_csv:
    resources:
        mem_mb=10000
    input:
        cf_actionables_csv
    output:
        vcf = 'actionables.vcf',
        bed = 'actionables.bed',
        dummy = 'dummy.txt'          # this dummy output is needed to bind this rule with the rule 'make_pileup' that needs to wait for the output of this rule
    shell:
        """
        bash {cf_snakedir}/scripts/csv_to_bed.sh "{input}" "{output.bed}"
        python3 {cf_snakedir}/scripts/csv_to_vcf.py "{input}" "{output.vcf}"
        echo "Hello there" > "{output.dummy}"
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
        bamlist = cf_bamlist_file,
        dummy = 'dummy.txt'
    output:
        'pileup/pileup.txt'
    params:
        bed = 'actionables.bed',
        rule_reference_genome = cf_reference_genome
    shell:
        """
        bash {cf_snakedir}/scripts/pileup.sh "{input.bamlist}" "{params.bed}" "{params.rule_reference_genome}" "{output}"
        """


rule process_pileup:
    resources:
        mem_mb = 100000
    input:
        actionables = cf_actionables_csv,
        pileup = 'pileup/pileup.txt',
        bamlist_file = cf_bamlist_file
    output:
        'processed_pileup/pileup.processed.csv'
    shell:
        """
        python3 {cf_snakedir}/scripts/process_pileup.py "{input.actionables}" "{input.pileup}" "{output}" "{input.bamlist_file}"
        """

