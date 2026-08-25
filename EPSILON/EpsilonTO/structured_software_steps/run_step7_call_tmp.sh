#!/usr/bin/bash

data_call="calling_data_dir"
error_model="trained_model.weights.h5"

python3 scripts/call_tmp.py --calling_data_dir $data_call --error_model $error_model



