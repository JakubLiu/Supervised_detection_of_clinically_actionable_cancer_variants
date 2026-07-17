#!/usr/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


"$SCRIPT_DIR/scripts/Epsilon_Fit.sh" \
    --input NEGATIVE_CONTROL_DATA3.txt \
    --noise_level 0.01 \
    --alt_mode generic \
    --output_model alt_generic_model