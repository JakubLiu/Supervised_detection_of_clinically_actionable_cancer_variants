import numpy as np
import tensorflow as tf
import random
from poligon import *  # your Epsilon class and submodels

# ----------------------------
# Simulation parameters
# ----------------------------
num_samples = 100               # total number of samples
features_per_read = 5
max_reads_possible = 500  # just the maximum number reads that can be simulated
padded_size = 10000     # the padding size to which the data will be padded
window_size = 21
num_unique_bases = 5
num_genomic_features = 10
num_features_per_sample = 10
base_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

# ----------------------------
# Function to pad reads and create masks
# ----------------------------

# this function takes as input a raw dataset and returns the padded dataset and the padding mask
def pad_reads_and_create_mask(reads, padded_size):
    coverage = reads.shape[0]
    features = reads.shape[1]
    padded = np.zeros((padded_size, features), dtype=np.float32)
    padded[:coverage, :] = reads
    mask = np.zeros((padded_size,), dtype=np.float32)
    mask[:coverage] = 1.0
    return padded, mask

# ----------------------------
# Simulate dataset
# ----------------------------


"""
- loop over the number of samples to simulate
- for each sample return the padded dataset and the corresponding padded mask
- concat them inot a big matrix of data and a big matrix of padding masks
- offcourse also simulate the data that does not require padding 
"""
x_free_sample_padded = []
x_free_sample_mask = []
x_free_genomic = []
x_instructive_genomic = []
x_instructive_sample = []
n_true = []  # just coverage put into another tensor
k_true = []  # the alt read count (I allow for at most 5% of the reads to be alternative)

for _ in range(num_samples):
    # Variable-length reads
    coverage = random.randint(50, max_reads_possible)
    reads = np.random.randn(coverage, features_per_read).astype(np.float32)
    padded, mask = pad_reads_and_create_mask(reads, padded_size)
    x_free_sample_padded.append(padded)
    x_free_sample_mask.append(mask)
    
    # Genomic sequence
    seq = ''.join(random.choice(['A','T','G','C','N']) for _ in range(window_size))
    x_free_genomic.append([base_to_int[b] for b in seq])
    
    # Instructive features
    x_instructive_genomic.append(np.random.randn(num_genomic_features).astype(np.float32))
    x_instructive_sample.append(np.random.randn(num_features_per_sample).astype(np.float32))
    
    # Target
    n_true.append(coverage)
    k_true.append(random.randint(0, int(coverage * 0.05)))


"""
Here the full dataset is created, it will get batched later.
"""
# Convert lists to arrays
x_free_sample_padded = np.stack(x_free_sample_padded, axis=0)
x_free_sample_mask = np.stack(x_free_sample_mask, axis=0)
x_free_genomic = np.array(x_free_genomic, dtype=np.int32)
x_instructive_genomic = np.stack(x_instructive_genomic, axis=0)
x_instructive_sample = np.stack(x_instructive_sample, axis=0)
n_true = np.stack(n_true, axis=0)
k_true = np.stack(k_true, axis=0)
y_true = y_targets = np.stack([n_true, k_true], axis=1)

# ----------------------------
# Instantiate Epsilon model
# ----------------------------
epsilon_model = Epsilon(
    intructive_genomic_input_dim=num_genomic_features,
    instructive_genomic_hidden_dim1=64,
    instructive_genomic_hidden_dim2=32,
    instructive_genomic_output_dim=16,

    instructive_sample_input_dim=num_features_per_sample,
    instructive_sample_hidden_dim1=64,
    instructive_sample_hidden_dim2=32,
    instructive_sample_output_dim=16,

    free_sample_hidden_dim1_psi=64,
    free_sample_hidden_dim2_psi=32,
    free_sample_output_dim_psi=16,
    free_sample_hidden_dim1_phi=64,
    free_sample_hidden_dim2_phi=32,
    free_sample_output_dim_phi=16,

    genomic_context_size=window_size,
    num_unique_bases=num_unique_bases,
    base_embedding_dim=128,
    dense_dim=32,
    num_heads=4,

    final_proj_dim1=32,
    final_proj_dim2=16
)

# ----------------------------
# Wrap in Functional API for dictionary input
# ----------------------------
instructive_genomic_input = tf.keras.Input(shape=(num_genomic_features,), name="instructive_genomic_input")
instructive_sample_input = tf.keras.Input(shape=(num_features_per_sample,), name="instructive_sample_input")
free_sample_padded_input = tf.keras.Input(shape=(padded_size, features_per_read), name="free_sample_padded_input")
free_genomic_input = tf.keras.Input(shape=(window_size,), name="free_genomic_input")
free_sample_padding_mask = tf.keras.Input(shape=(padded_size,), name="free_sample_padding_mask")

output = epsilon_model(
    instructive_genomic_input,
    instructive_sample_input,
    free_sample_padded_input,
    free_genomic_input,
    free_sample_padding_mask
)

functional_model = tf.keras.Model(
    inputs={
        'instructive_genomic_input': instructive_genomic_input,
        'instructive_sample_input': instructive_sample_input,
        'free_sample_padded_input': free_sample_padded_input,
        'free_genomic_input': free_genomic_input,
        'free_sample_padding_mask': free_sample_padding_mask
    },
    outputs=output
)

functional_model.compile(
    optimizer=tf.keras.optimizers.Adam(1e-3),
    loss=BinomialNegLogLoss(pseudocount=1e-8, normalize_by_coverage=False),
    metrics=['mae']
)

# ----------------------------
# Create tf.data.Dataset for batching
# ----------------------------


"""
Here the full dataset is getting batched.
"""

batch_size = 8
dataset = tf.data.Dataset.from_tensor_slices((
    {
        'instructive_genomic_input': x_instructive_genomic,
        'instructive_sample_input': x_instructive_sample,
        'free_sample_padded_input': x_free_sample_padded,
        'free_genomic_input': x_free_genomic,
        'free_sample_padding_mask': x_free_sample_mask
    },
    y_true
))

dataset = dataset.shuffle(buffer_size=num_samples).batch(batch_size)

# ----------------------------
# Train the model
# ----------------------------
functional_model.fit(dataset, epochs=200)

print("Training finished!")