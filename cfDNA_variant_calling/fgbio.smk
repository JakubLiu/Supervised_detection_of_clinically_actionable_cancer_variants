import numpy as np  


R1_paths = list(
    np.loadtxt('R1.paths.sorted.txt', dtype = str)
)

R2_paths = list(
    np.loadtxt('R2.paths.sorted.txt', dtype = str)
)

R1_names = list(
    np.loadtxt('R1.names.sorted.txt', dtype = str)
)

R2_names = list(
    np.loadtxt('R2.names.sorted.txt', dtype = str)
)


R1_R2_dict = {}
R1_name_R1_path_dict = {}

for i in range(0, len(R1_paths)):
    R1_R2_dict[R1_paths[i]] = R2_paths[i]
    R1_name_R1_path_dict[R1_names[i]] = R1_paths[i]



rule all:
    input:
        expand('output/{sample}.consensus.mapped.unfiltered.bam', sample = R1_names),
        expand('output/{sample}.consensus.mapped.filtered.bam', sample = R1_names)


rule fastq_to_ubam:
    resources:
        mem_mb = 32000
    input:
        R1 = lambda wc: R1_name_R1_path_dict[wc.sample],
        R2 = lambda wc: R1_R2_dict[R1_name_R1_path_dict[wc.sample]]
    output:
        'output/{sample}.unmapped.bam'
    params:
        read_struct = '3M2S+T 3M2S+T',
        samplename = lambda wc: wc.sample
    shell:
        """
        fgbio -Xmx1g --compression 1 --async-io FastqToBam \
            --input {input.R1} {input.R2} \
            --read-structures {params.read_struct} \
            --sample {params.samplename} \
            --library library1 \
            --output {output}
        """


rule ubam_to_mapped_bam:
    resources:
        mem_mb = 32000
    input:
        'output/{sample}.unmapped.bam'
    output:
        'output/{sample}.mapped.bam'
    params:
        ref = 'hs37d5.fa'
    shell:
        """
        samtools fastq {input} \
            | bwa mem -t 16 -p -K 150000000 -Y {params.ref} - \
            | fgbio -Xmx4g --compression 1 --async-io ZipperBams \
                --unmapped {input} \
                --ref {params.ref} \
                --output {output}
        """


rule sort_mapped_bam:
    resources:
        mem_mb = 32000
    input:
        'output/{sample}.mapped.bam'
    output:
        'output/{sample}.mapped.sorted.bam'
    threads:
        8
    shell:
        """
        samtools sort --template-coordinate --threads {threads} \
            -o {output} \
            {input}
        """


rule group_bam:
    resources:
        mem_mb = 32000
    input:
        'output/{sample}.mapped.sorted.bam'
    output:
        bam = 'output/{sample}.mapped.sorted.grouped.bam',
        hist = 'output/{sample}.hist.txt'
    shell:
        """
        fgbio -Xmx8g --compression 1 --async-io GroupReadsByUmi \
            --input {input} \
            --strategy paired \
            --edits 1 \
            --output {output.bam} \
            --family-size-histogram {output.hist}
        """


rule consensus_bam:
    resources:
        mem_mb = 32000
    input:
        'output/{sample}.mapped.sorted.grouped.bam'
    output:
        'output/{sample}.consensus.unmapped.bam'
    params:
        ref = 'hs37d5.fa'
    threads:
        8
    shell:
        """
        fgbio -Xmx4g --compression 0 CallDuplexConsensusReads \
            --input {input} \
            --output /dev/stdout \
            --min-reads 3 \
            --min-input-base-quality 20 \
            --threads {threads} \
            |  fgbio -Xmx8g --compression 1 FilterConsensusReads \
                --input /dev/stdin \
                --output {output} \
                --ref {params.ref} \
                --min-reads 3 \
                --min-base-quality 45 \
                --max-base-error-rate 0.2
        """



rule map_consensus_bam:
    resources:
        mem_mb = 32000
    input:
        'output/{sample}.consensus.unmapped.bam'
    output:
        unfiltered = 'output/{sample}.consensus.mapped.unfiltered.bam',
        filtered = 'output/{sample}.consensus.mapped.filtered.bam'
    params:
        ref = 'hs37d5.fa'
    threads:
        8
    shell:
        """
        # Re-align the consensus reads
        samtools fastq {input} \
        | bwa mem -t {threads} -p -K 150000000 -Y {params.ref} - \
        | fgbio -Xmx4g --compression 1 --async-io ZipperBams \
            --unmapped {input} \
            --ref {params.ref} \
            --tags-to-reverse Consensus \
            --tags-to-revcomp Consensus \
            --output {output.unfiltered}

        # Filter and sort the consensus reads
        fgbio -Xmx8g --compression 0 FilterConsensusReads \
        --input {output.unfiltered} \
        --output /dev/stdout \
        --ref {params.ref} \
        --min-reads 3 \
        --min-base-quality 45 \
        --max-base-error-rate 0.2 \
        | samtools sort --threads 8 -o {output.filtered} --write-index
        """