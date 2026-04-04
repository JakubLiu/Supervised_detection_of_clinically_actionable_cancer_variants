import numpy as np
import pandas as pd
import tensorflow as tf
import random
from poligon import *  # for Epsilon

def pad_reads_and_create_mask(reads, padded_size):
    coverage = reads.shape[0]
    features = reads.shape[1]
    padded = np.zeros((padded_size, features), dtype=np.float32)
    padded[:coverage, :] = reads
    mask = np.zeros((padded_size,), dtype=np.float32)
    mask[:coverage] = 1.0
    return padded, mask

base_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
data_samples = pd.read_csv('../make_training_data/data_negative_controls.csv')
data_reads = pd.read_csv('../make_training_data/per_read_data.csv')


data_reads['chrom'] = data_reads['chrom'].astype(str)
data_samples['pos'] = data_samples['pos'].astype(int)
data_reads["strand"] = data_reads["strand"].map({"+": 1, "-": 0})

padded_size = 10000

x_free_sample_padded = []
x_free_sample_mask = []
x_free_genomic = []
x_instructive_genomic = []
x_instructive_sample = []
n_true = []  # just coverage put into another tensor
k_true = []  # the alt read count (I allow for at most 5% of the reads to be alternative)


for i in range(0, data_samples.shape[0]):
    
    print(str(i/data_samples.shape[0]*100) + '%', flush = True)
    current_chrom = data_samples.iloc[i,0]
    current_pos = data_samples.iloc[i,1]
    current_sampleID = data_samples.iloc[i,7]
    
    # extract the elements of the per read data that correspond to the given locus and the given patient
    corresponding_per_read_data = data_reads[(data_reads['chrom'] == current_chrom) &
                                             (data_reads['pos'] == current_pos) &
                                             (data_reads['sampleID'] == current_sampleID)]
    
    # padd the data and compute the padding mask
    padded, mask = pad_reads_and_create_mask(corresponding_per_read_data.iloc[:,3:], padded_size)
    x_free_sample_padded.append(padded)
    x_free_sample_mask.append(mask)

    raw_genomic_context = list(data_samples.iloc[i,-1])[(5000-10):(5000+10+1)]  # extract only a 21bp window around the locus of interest
    x_free_genomic.append([base_to_int[x] for x in raw_genomic_context])
    x_instructive_genomic.append([float(x) for x in list(data_samples.loc[i, ["GGC_upstream","CGG_downstream","GGT_upstream","TGG_downstream","GC_content","homo_percentage"]])])  # remove some homo feature due to almonst no variance
    x_instructive_sample.append([float(x) for x in list(data_samples.loc[i, ["median_mapQ","median_baseQ", "mismatches"]])])
    k_true.append(float(data_samples.loc[i,'specific_alt_counts']))  # alt counts
    n_true.append(float(data_samples.loc[i, 'general_alt_counts']) + float(data_samples.loc[i,'ref_counts']))  # coverage
  


x_free_sample_padded = np.stack(x_free_sample_padded, axis=0)
x_free_sample_mask = np.stack(x_free_sample_mask, axis=0)
x_free_genomic = np.array(x_free_genomic, dtype=np.int32)
x_instructive_genomic = np.stack(x_instructive_genomic, axis=0)
x_instructive_sample = np.stack(x_instructive_sample, axis=0)
n_true = np.stack(n_true, axis=0)
k_true = np.stack(k_true, axis=0)
y_true = np.stack([n_true, k_true], axis=1)


train_mask = np.random.uniform(0,1, y_true.shape[0]) < 0.9

# Free sample (padded + mask)
x_free_sample_padded_train = x_free_sample_padded[train_mask]
x_free_sample_padded_val   = x_free_sample_padded[~train_mask]

x_free_sample_mask_train = x_free_sample_mask[train_mask]
x_free_sample_mask_val   = x_free_sample_mask[~train_mask]

# Free genomic
x_free_genomic_train = x_free_genomic[train_mask]
x_free_genomic_val   = x_free_genomic[~train_mask]

# Instructive
x_instructive_genomic_train = x_instructive_genomic[train_mask]
x_instructive_genomic_val   = x_instructive_genomic[~train_mask]

x_instructive_sample_train = x_instructive_sample[train_mask]
x_instructive_sample_val   = x_instructive_sample[~train_mask]

# Targets
y_true_train = y_true[train_mask]
y_true_val   = y_true[~train_mask]


np.savez("training_data.npz",
    x_instructive_genomic_train=x_instructive_genomic_train,
    x_instructive_sample_train=x_instructive_sample_train,
    x_free_sample_padded_train=x_free_sample_padded_train,
    x_free_genomic_train=x_free_genomic_train,
    x_free_sample_mask_train=x_free_sample_mask_train,
    y_true_train=y_true_train
)



np.savez("validation_data.npz",
    x_instructive_genomic_val=x_instructive_genomic_val,
    x_instructive_sample_val=x_instructive_sample_val,
    x_free_sample_padded_val=x_free_sample_padded_val,
    x_free_genomic_val=x_free_genomic_val,
    x_free_sample_mask_val=x_free_sample_mask_val,
    y_true_val=y_true_val
)

print('done.')
