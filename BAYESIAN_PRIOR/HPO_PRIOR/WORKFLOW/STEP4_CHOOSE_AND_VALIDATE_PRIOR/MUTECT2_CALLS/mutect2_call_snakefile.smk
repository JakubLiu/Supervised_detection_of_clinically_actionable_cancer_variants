

normal_bams = 'normals.txt'
tumor_bams = 'tumors.txt'
samples = 'samplenames.txt'



sample_to_tumor_dict = {}   # BIH_166 -> path
sample_to_normal_dict = {}  # BIH_166 -> path
sample_to_tumor_name_dict = {}  # BIH_166 -> BIH_166_T1
sample_to_normal_name_dict = {}  # BIH_166 -> BIH_166_N1

with open(samples, 'r') as s, open(normal_bams, 'r') as n, open(tumor_bams, 'r') as t:

    for samplename, normal, tumor in zip(s,n,t):
        samplename = samplename.strip()
        normal = normal.strip()
        tumor = tumor.strip()

        sample_to_normal_dict[samplename] = normal
        sample_to_tumor_dict[samplename]  = tumor
        sample_to_normal_name_dict[samplename] = normal.split('/')[-1].split('.')[1]
        sample_to_tumor_name_dict[samplename] = tumor.split('/')[-1].split('.')[1]



print(sample_to_normal_name_dict)
print('-'*50)
print(sample_to_tumor_name_dict)
print('-'*50)
print(sample_to_normal_dict)
print('-'*50)
print(sample_to_tumor_dict)


rule all:
    input:
        expand('mutect2_PASS/{sample}.PASS.vcf.gz', sample = list(sample_to_tumor_dict.keys()))


rule mutect2_call:
    resources:
        mem_mb = 100000
    input:
        tumor = lambda wc: sample_to_tumor_dict[wc.sample],
        normal = lambda wc: sample_to_normal_dict[wc.sample]
    output:
        'mutect2_raw_call/{sample}.vcf.gz'
    params:
        reference_genome = '/data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/TUNE_PRIOR/STEP2_SEARCH_FOR_CASES/hs37d5.fa',
        normal_sample_name = lambda wc: sample_to_normal_name_dict[wc.sample]
    shell:
        """
        bash /data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/TUNE_PRIOR/STEP5_VALIDATE/MUTECT2_CALLS/scripts/mutect2_call.sh \
                                                                                                                        "{input.tumor}" \
                                                                                                                        "{input.normal}" \
                                                                                                                        "{params.normal_sample_name}" \
                                                                                                                        "{params.reference_genome}" \
                                                                                                                        "{output}"
        """


rule mutect2_filter:
    resources:
        mem_mb=100000
    input:
        'mutect2_raw_call/{sample}.vcf.gz'
    output:
        'mutect2_flagged/{sample}.flagged.vcf.gz'
    params:
        reference_genome = '/data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/TUNE_PRIOR/STEP2_SEARCH_FOR_CASES/hs37d5.fa'
    shell:
        """
        bash /data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/TUNE_PRIOR/STEP5_VALIDATE/MUTECT2_CALLS/scripts/mutect2_filter.sh \
                                                                                                                            "{params.reference_genome}" \
                                                                                                                            "{input}" \
                                                                                                                            "{output}"
        """


rule keep_PASS:
    resources:
        mem_mb=100000
    input:
        'mutect2_flagged/{sample}.flagged.vcf.gz'
    output:
        'mutect2_PASS/{sample}.PASS.vcf.gz'
    shell:
        """
        bash /data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/TUNE_PRIOR/STEP5_VALIDATE/MUTECT2_CALLS/scripts/keep_PASS.sh \
                                                                                                                            "{input}" \
                                                                                                                            "{output}"
        """