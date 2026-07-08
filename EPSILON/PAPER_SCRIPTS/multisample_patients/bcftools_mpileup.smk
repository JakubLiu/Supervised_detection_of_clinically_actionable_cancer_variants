import numpy as np  

tumor_paths = list(
    np.loadtxt('paths_tumors.txt', dtype = str)
)

tumor_names = list(
    np.loadtxt('tumor_names.txt', dtype = str)
)

name_path_dict = {}
for i in range(0, len(tumor_paths)):
    name_path_dict[tumor_names[i]] = tumor_paths[i]


print(name_path_dict)

rule all:
    input:
        expand('bcftools_mpileup/{sample}.txt', sample = tumor_names)

rule bcftools_mpileup:
    resources:
        mem_mb = 32000
    input:
        lambda wc: name_path_dict[wc.sample]
    output:
        'bcftools_mpileup/{sample}.txt'
    params:
        ref = 'hs37d5.fa'
    shell:
        """
        # keep only sited that have >= 2 non-reference alleles

        bcftools mpileup -Ou -f {params.ref} {input} | \
        bcftools call -mv -Ou | \
        bcftools view -i 'FORMAT/AD[0:1] >= 2' -Ou | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' > {output}
        """
