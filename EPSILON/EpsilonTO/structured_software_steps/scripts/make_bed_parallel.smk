"""
How to run:

    snakemake --snakefile make_bed_parallel.smk --configfile make_bed_parallel_config.yaml --profile <str> --jobs <int>

"""




import os


 # #######################################################################################################################

refgen = config['reference_genome']
fai = config['reference_genome_fai']
bam = config['bamfile']
chromosome = config['chromosome']
chunk_size = config['chunk_size']
minimum_coverage = config['min_coverage']
output_file = config['output']

# ###############################################################################################################################

if not os.path.exists(fai):
    raise FileNotFoundError(
        f"{fai} does not exist. Run: samtools faidx {refgen}"
    )

chrom_length = None

with open(fai) as f:
    for line in f:
        fields = line.rstrip().split("\t")
        if fields[0] == chromosome:
            chrom_length = int(fields[1])
            break

if chrom_length is None:
    raise ValueError(f"Chromosome {chromosome} not found in {fai}")


chunks = {}

start = 1  # samtools -r uses 1-based coordinates

while start <= chrom_length:
    end = min(start + chunk_size - 1, chrom_length)

    # filename-safe name
    chunk_name = f"{chromosome}_{start}_{end}"

    # samtools region
    region = f"{chromosome}:{start}-{end}"

    chunks[chunk_name] = region

    start = end + 1


# #############################################################################################################################################

rule all:
    input:
        output_file


rule samtools_mpileup_chunk:
    input:
        bam=bam
    output:
        "chunks/{chunk}.mpileup.txt"
    params:
        ref=refgen,
        region=lambda wildcards: chunks[wildcards.chunk]
    resources:
        mem_mb=32000
    shell:
        """
        mkdir -p chunks

        samtools mpileup \
            -f {params.ref} \
            -r {params.region} \
            {input.bam} \
            -aa \
            > {output}
        """


rule merge_mpileup:
    input:
        expand(
            "chunks/{chunk}.mpileup.txt",
            chunk=chunks.keys()
        )
    output:
        "merged.mpileup.txt"
    shell:
        """
        cat {input} > {output}
        """

rule process_pileup:
    resources:
        mem_mb = 32000
    input:
        "merged.mpileup.txt"
    output:
        output_file
    params:
        min_coverage = minimum_coverage
    shell:
        """
        awk -F '\t' '$4 > {params.min_coverage} {{print}}' {input} | awk -F '\t' '{{print $1, $2-1, $2, $3}}' | sed 's/ /\t/g' > {output}
        """