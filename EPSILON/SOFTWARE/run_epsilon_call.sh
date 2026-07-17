#!/usr/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# for bayes generic____________________________________________________________________________
#"$SCRIPT_DIR/scripts/Epsilon_call.sh" \
#    --model fitted_error_model/alt_generic_model.rds \
#    --input tumor_data/tumor_data3.txt \
#    --output generic.bayes.vcf \
#    --mode bayes_posterior \
#    --prior 0.3 \
#    --posterior_cutoff 0.5 \
#   --alt_mode generic


# for binomial test generic _________________________________________________________________
"$SCRIPT_DIR/scripts/Epsilon_call.sh" \
    --model fitted_error_model/alt_generic_model.rds \
    --input tumor_data/tumor_data3.txt \
    --output generic.binomial.vcf \
    --mode binomial_test \
    --fdr_method BH \
    --alpha 0.05 \
   --alt_mode generic



# for bayes specific____________________________________________________________________________
#"$SCRIPT_DIR/scripts/Epsilon_call.sh" \
#    --model fitted_error_model/alt_specific_model.rds \
#    --input tumor_data/tumor_data2.txt \
#    --output specific.bayes.vcf \
#    --mode bayes_posterior \
#    --prior 0.3 \
#    --posterior_cutoff 0.5 \
#   --alt_mode specific


# for binomial test specific _________________________________________________________________
#"$SCRIPT_DIR/scripts/Epsilon_call.sh" \
#    --model fitted_error_model/alt_specific_model.rds \
#    --input tumor_data/tumor_data2.txt \
#    --output specific.binomial.vcf \
#    --mode binomial_test \
#    --fdr_method BH \
#    --alpha 0.05 \
#    --alt_mode specific