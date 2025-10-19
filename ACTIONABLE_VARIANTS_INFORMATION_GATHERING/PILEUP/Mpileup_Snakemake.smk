# config file
reference_genome = config['reference_genome']
bamfile_list = config['bamfile_list']


rule all:
    input:
        'processed_mpileup/pileup_processed.txt'


rule samtools_mpileup:
    input:
        ref_gen = reference_genome,
        bamlist = bamfile_list
    output:
        'samtools_mpileup/pileup_raw.txt'

    resources:
        mem='60G',
        time='20:00:00'

    shell:
        '''
        samtools mpileup -f {input.ref_gen} -b {input.bamlist} -s > {output}
        '''

rule postprocess_mpileup:
    input:
        'samtools_mpileup/pileup_raw.txt'
    output:
        'processed_mpileup/pileup_processed.txt'

    resources:
        mem='60G',
        time='40:00:00'

    script:
        "scripts/postprocess_mpileup.py"
