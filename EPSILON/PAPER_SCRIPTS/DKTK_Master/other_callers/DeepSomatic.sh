#!/usr/bin/bash

# =================================== usage =================================
# ./DeepSomatic.sh \
#        /data/processed/ \ ---> dir to bam files
#        /data/precomputed/ \   --> dir to reference genome
#        T1-DNA1-WES1.sorted.bam \  --> tumor bam
#        N1-DNA1-WES1.sorted.bam\  --> normal bam
#        "T1-DNA1-WES1" \  --> tumor samplename
#        "N1-DNA1-WES1" \  --> normal samplename
#        test_output.vcf.gz\  --> output vcf
#        test_output.gvcf.gz \  --> output genomic vcf
#        Civic.small.bed   --> bed file
# ===========================================================================


BIN_VERSION="1.10.0"

# singularity pull docker://google/deepsomatic:"${BIN_VERSION}" --> run this if not done before

bamdir="$1"
refdir="$2"
tumor_bam="$3"
normal_bam="$4"
tumor_samplename="$5"
normal_samplename="$6"
output_vcf="$7"
output_gvcf="$8"
bed="$9"




singularity run \
-B /usr/lib/locale/:/usr/lib/locale/ \
-B $HOME:$HOME \
-B $refdir:/precomputed \
-B $bamdir:/bams \
docker://google/deepsomatic:"${BIN_VERSION}" \
run_deepsomatic \
--model_type=WES \
--ref=/precomputed/hs37d5.fa \
--reads_normal=/bams/$normal_bam \
--reads_tumor=/bams/$tumor_bam \
--output_vcf=$output_vcf \
--output_gvcf=$output_gvcf \
--sample_name_tumor=$tumor_samplename \
--sample_name_normal=$normal_samplename \
--num_shards=1 \
--logging_dir=logs \
--intermediate_results_dir intermediate_results_dir \
--regions="$(awk '{printf "%s:%d-%d ", $1, $2, $3} END {print ""}' $bed)"
