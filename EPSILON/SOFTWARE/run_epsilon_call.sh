#!/usr/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# for bayes____________________________________________________________________________
#"$SCRIPT_DIR/scripts/Epsilon_call.sh" \
#    --model alt_specific_fitted_model.rds \
#    --input tumor_data/tumor_data2.txt \
#    --output test.bayes.vcf \
#    --mode bayes_posterior \
#    --prior 0.3 \
#    --posterior_cutoff 0.5 \
#    --alt_mode specific


# for binomial test _________________________________________________________________
"$SCRIPT_DIR/scripts/Epsilon_call.sh" \
    --model fitted_error_model/alt_specific_fitted_model.rds \
    --input tumor_data/tumor_data2.txt \
    --output test.binomial.vcf \
    --mode binomial_test \
    --fdr_method BH \
    --alpha 0.05 \
    --alt_mode specific