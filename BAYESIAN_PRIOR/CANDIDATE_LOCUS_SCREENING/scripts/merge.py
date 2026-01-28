import pandas as pd
import numpy as np
import sys

tumor = sys.argv[1]  # processed tumor pileup
normal = sys.argv[2]   # processed normal pileup
output_merged = sys.argv[3]    # merged output df


tumor = pd.read_csv(tumor)
normal = pd.read_csv(normal)

merged = pd.merge(tumor, normal, on = ['chrom', 'pos', 'ref'], how = 'inner')

merged = merged.loc[merged['total_alt_count_x'] > 0,]
merged = merged.loc[merged['total_alt_count_y'] < 10,]    # this is just a very rough preliminary filter

merged.to_csv(output_merged, index = False)

