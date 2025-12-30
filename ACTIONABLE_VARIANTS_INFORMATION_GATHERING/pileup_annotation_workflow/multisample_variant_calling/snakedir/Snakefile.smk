import os
import re
 
# helper functions =======================================================================================================
def extract_basename(path):
    filename = os.path.basename(path)
    name = re.sub(r'\..*\.bam$', '', filename)
    return name
 
# ======================================================================================================================
 
 
# get data from the config file ======================================================================================
cf_snakedir = config['snakedir']
cf_reference_genome = config['reference_genome']
cf_snpEff_dir = config['snpEff_dir']
cf_snpEff_annotation_file = config['snpEff_annotation_file']
# =================================================================================================================


rule all:
    input:
      'pileup/pileup.processed.csv',
      'annotation/annotated.vcf'


rule make_vcf_to_annotate:
    resources:
        mem_mb=config['resources_mb']['make_vcf_to_annotate']
    input:
        'init_input_files/multi_vcf_compilation.txt'
    output:
        loci = 'variants_loci/loci.txt',
        vcf = 'annotation/vcf_to_annotate.vcf'
    shell:
        """
        bash {cf_snakedir}/scripts/extract_loci.sh "{input}" tmp.txt
        python3 {cf_snakedir}/scripts/make_vcf.py "{input}" "{output.vcf}"
        sed 's/;/,/g' tmp.txt |  tail -n +2 > "{output.loci}"
        rm tmp.txt
        """
 


rule annotate:
    resources:
        mem_mb=config['resources_mb']['annotate']
    input:
        'annotation/vcf_to_annotate.vcf'
    output:
        'annotation/annotated.vcf'
    params:
        rule_annotation_file = cf_snpEff_annotation_file,
        rule_snpEff_dir = cf_snpEff_dir
    shell:
        """
        bash {cf_snakedir}/scripts/annotate.sh "{input}" "{params.rule_annotation_file}" "{output}" "{params.rule_snpEff_dir}"
        """



rule pileup:
    resources:
        mem_mb=config['resources_mb']['pileup']
    input:
        bamlist = config['bamlist_file']
    output:
        'pileup/raw_pileup.txt'
    params:
        rule_reference_genome = cf_reference_genome
    shell:
        """
        awk -F ';' -v OFS='\t' '{{print $1, $2-1, $2}}' {cf_snakedir}/init_input_files/multi_vcf_compilation.txt | tail -n +2 > tmp2.txt
        bash {cf_snakedir}/scripts/pileup.sh "{input.bamlist}" tmp2.txt "{params.rule_reference_genome}" "{output}"
        rm tmp2.txt
        """




rule process_pileup:
    resources:
        mem_mb = config['resources_mb']['process_pileup']
    input:
        loci = 'variants_loci/loci.txt',
        pileup = 'pileup/raw_pileup.txt',
        bamlist_file = config['bamlist_file'],
        multi_vcf = 'init_input_files/multi_vcf_compilation.txt'
    output:
        'pileup/pileup.processed.csv'
    shell:
        """
        python3 {cf_snakedir}/scripts/process_pileup.py "{input.loci}" "{input.pileup}" "{output}" \
                                                        "{input.bamlist_file}" "{input.multi_vcf}"
        """

