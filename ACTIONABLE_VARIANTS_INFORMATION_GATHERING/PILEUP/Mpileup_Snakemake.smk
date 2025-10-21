# config file
reference_genome = config['reference_genome']
bamfile_list = config['bamfile_list']
regions_bed = config['regions_bed']


rule all:
    input:
        'processed_mpileup/pileup_processed_regions.txt'


rule samtools_mpileup:
    input:
        ref_gen = reference_genome,
        bamlist = bamfile_list,
        bedfile = regions_bed
    output:
        'samtools_mpileup/pileup_raw_regions.txt'

    resources:
        mem='60G',
        time='20:00:00'

    shell:
        '''
        samtools mpileup -f {input.ref_gen} -b {input.bamlist} -l {input.bedfile} -s > {output}
        '''


rule postprocess_mpileup:
    input:
        'samtools_mpileup/pileup_raw_regions.txt'
    output:
        'processed_mpileup/pileup_processed_regions.txt'

    resources:
        mem='60G',
        time='40:00:00'

    script:
        "scripts/postprocess_mpileup_optimized.py"
