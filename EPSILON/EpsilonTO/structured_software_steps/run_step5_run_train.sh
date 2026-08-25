#!/usr/bin/bash







training_data_dir="npz_output_test"

python3 scripts/train.py \
    --data_dir "$training_data_dir" \
    --max_epochs 1000 \
    --zero_keep_prop_multiplier 10.0 \
    --optimizer adam \
    --learning_rate 1e-4 \
    --patience 3 \
    --max_error_rate 0.05 \
    --output_weights fitted.weights.h5



    
