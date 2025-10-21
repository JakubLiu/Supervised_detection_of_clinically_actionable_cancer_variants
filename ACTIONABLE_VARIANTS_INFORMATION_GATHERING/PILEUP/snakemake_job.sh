#!/usr/bin/bash

snakemake --snakefile Mpileup_Snakemake.smk --configfile config.yaml --profile=cubi-v1 --jobs 1
