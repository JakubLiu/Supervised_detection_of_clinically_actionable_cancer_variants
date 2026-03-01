import os

BAMS = []  # just the sample names, not the paths

with open('bams_to_mutate_train.txt', 'r') as f:
    for line in f:
        BAMS.append(os.path.basename(line).replace('\n', '').replace('.bam', ''))



bam_dict = {}
with open('bams_to_mutate_train.txt') as f:
    for line in f:
        path = line.strip()
        sample_name = os.path.basename(path).replace('.bam','')
        bam_dict[sample_name] = path




rule all:
    input:
        sorted_ = expand('inserted_SNV/{sample}.sorted.bam', sample = BAMS),
        tmp = expand('inserted_SNV/{sample}.mutated.bam', sample = BAMS),
        final = expand('inserted_SNV/{sample}.sorted.bam', sample = BAMS)


rule insert_SNV:
    resources:
        mem_mb=100000
    input:
        bam=lambda wildcards: bam_dict[wildcards.sample]
    output:
        sorted_ = temp('inserted_SNV/{sample}.sorted.bam'),
        tmp = temp("inserted_SNV/{sample}.mutated.bam"),
        final = "inserted_SNV/{sample}.mutated.sorted.bam"
    params:
        reference_genome='/hs37d5.fa',
        varfile='SNV_location_file.txt'
    shell:
        """
        samtools sort "{input.bam}" -o "{output.sorted_}"
        samtools index "{output.sorted_}"

        python3 /addsnv.py \
                           -v "{params.varfile}" -r "{params.reference_genome}" \
                           -f "{output.sorted_}" -o "{output.tmp}" \
                           --mindepth 1 --force \
        
        samtools sort "{output.tmp}" -o "{output.final}"
        samtools index "{output.final}"
        """
