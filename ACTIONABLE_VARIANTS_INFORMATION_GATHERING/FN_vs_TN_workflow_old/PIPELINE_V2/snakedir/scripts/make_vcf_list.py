#!/usr/bin/python

import sys
import os

output_file = sys.argv[1]
vcf_files = sys.argv[2:]

with open(output_file, "w") as f:
    for vcf in vcf_files:
        f.write(os.path.abspath(vcf) + "\n")