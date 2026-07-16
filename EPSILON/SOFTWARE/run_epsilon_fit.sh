#!/usr/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


"$SCRIPT_DIR/scripts/Epsilon_Fit.sh" \
    --input negative_control_data/NEGATIVE_CONTROL_DATA.txt \
    --noise_level 0.01 \
    --alt_mode specific \
    --output_model alt_specific_fitted_model2