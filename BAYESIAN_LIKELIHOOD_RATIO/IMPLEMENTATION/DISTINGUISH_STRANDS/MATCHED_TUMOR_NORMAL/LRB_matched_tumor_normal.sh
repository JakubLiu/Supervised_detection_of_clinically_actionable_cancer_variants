#!/usr/bin/bash

PADDING_UPSTREAM=10
PADDING_DOWNSTREAM=10
PSEUDOCOUNT=0.00001
PRIOR=0.005
POSTERIOR_CUTOFF=0.5
NORMAL_POSTERIOR_EVIDENCE_WARNING_THRESHOLD=0.2

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in

    --tumor_bam)
      TUMOR_BAM="$2"
      shift 2
      ;;

    --matched_normal_bam)
      MATCHED_NORMAL_BAM="$2"
      shift 2
      ;;

    --negative_control_bamlist)
      NEGATIVE_CONTROL_BAMLIST="$2"
      shift 2
      ;;

    --chromosome)
      CHROMOSOME="$2"
      shift 2
      ;;

    --start)
      START="$2"
      shift 2
      ;;

    --stop)
      STOP="$2"
      shift 2
      ;;
    
    --ref_allele)
      REF_ALLELE="$2"
      shift 2
      ;;

    --alt_allele)
      ALT_ALLELE="$2"
      shift 2
      ;;

    --prior)
      PRIOR="$2"
      shift 2
      ;;

    --posterior_cutoff)
      POSTERIOR_CUTOFF="$2"
      shift 2
      ;;

    --pseudocount)
      PSEUDOCOUNT="$2"
      shift 2
      ;;

    --output_call_file)
      OUTPUT_CALL_FILE="$2"
      shift 2
      ;;

    --min_mapQ)
      MIN_MAPQ="$2"
      shift 2
      ;;

    --min_baseQ)
      MIN_BASEQ="$2"
      shift 2
      ;;

    --output_read_annotation_file)
      OUTPUT_READ_ANNOTATION_FILE="$2"
      shift 2
      ;;

    --reference_genome)
      REFERENCE_GENOME="$2"
      shift 2
      ;;

    --padding_upstream)
      PADDING_UPSTREAM="$2"
      shift 2
      ;;

    --padding_downstream)
      PADDING_DOWNSTREAM="$2"
      shift 2
      ;;

    --output_genomic_context_file)
      OUTPUT_GENOMIC_CONTEXT_FILE="$2"
      shift 2
      ;;

    --normal_posterior_evidence_threshold)
      NORMAL_POSTERIOR_EVIDENCE_WARNING_THRESHOLD="$2"
      shift 2
      ;;

  esac
done





Rscript LRB_matched_tumor_normal.R \
                    "$TUMOR_BAM" \
                    "$MATCHED_NORMAL_BAM" \
                    "$NEGATIVE_CONTROL_BAMLIST" \
                    "$CHROMOSOME" \
                    "$START" \
                    "$STOP" \
                    "$REF_ALLELE" \
                    "$ALT_ALLELE" \
                    "$PRIOR" \
                    "$POSTERIOR_CUTOFF" \
                    "$PSEUDOCOUNT" \
                    "$OUTPUT_CALL_FILE" \
                    "$MIN_MAPQ" \
                    "$MIN_BASEQ" \
                    "$OUTPUT_READ_ANNOTATION_FILE" \
                    "$REFERENCE_GENOME" \
                    "$PADDING_UPSTREAM" \
                    "$PADDING_DOWNSTREAM" \
                    "$OUTPUT_GENOMIC_CONTEXT_FILE" \
                    "$NORMAL_POSTERIOR_EVIDENCE_WARNING_THRESHOLD"
