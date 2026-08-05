#!/usr/bin/bash







DATA_DIR="npz_output_test"

python3 train.py \
    --data_dir "$DATA_DIR" \
    --max_epochs 100 \
    --zero_keep_prop_multiplier 10.0 \
    --optimizer adam \
    --learning_rate 1e-4 \
    --patience 10 \
    --max_error_rate 0.05