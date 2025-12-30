import pandas as pd
import sys

processed_pileup = sys.argv[1]
multi_vcf = sys.argv[2]
output = sys.argv[3]
 

processed_pileup = pd.read_csv(processed_pileup)
multi_vcf = pd.read_csv(multi_vcf, delimiter = ',', header = 0)

wide = processed_pileup.pivot(
    index=["chrom", "pos", "ref", "variant_alt"],
    columns="sample",
    values=[
        "coverage",
        "n_alt_reads",
        "n_variant_supporting_alt_reads",
        "translated_sequence"
    ]
)
 
 
wide.columns = [
    f"{metric}__{sample.split('.')[0]}"
    for metric, sample in wide.columns
]
 
wide = wide.reset_index()
 
# add a chrom;pos;ref;alt column to merge with multi_vcf
wide["variant"] = (
    wide["chrom"].astype(str) + ";" +
    wide["pos"].astype(str) + ";" +
    wide["ref"].astype(str) + ";" +
    wide["variant_alt"].astype(str)
)


merged = pd.merge(multi_vcf, wide, on="variant", how="inner")
 
merged = merged.drop(columns = ['variant'])
 

front = ['chrom', 'pos', 'ref', 'variant_alt']
cols = front + [c for c in merged.columns if c not in front]
merged = merged[cols]

 
merged.to_csv(output, index = False)
